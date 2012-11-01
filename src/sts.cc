#include "sts/likelihood.hpp"
#include "sts/moves.hpp"
#include "sts/particle.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <fstream>
#include <stack>
#include <memory>
#include <unordered_map>
#include <Bpp/Phyl/Model/GTR.h>
#include <Bpp/Phyl/Model/HKY85.h>
#include <Bpp/Phyl/Model/JCnuc.h>
#include <Bpp/Phyl/Model/JTT92.h>
#include <Bpp/Phyl/Model/TN93.h>
#include <Bpp/Phyl/Model/WAG01.h>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include "tclap/CmdLine.h"

#define _STRINGIFY(s) #s
#define STRINGIFY(s) _STRINGIFY(s)

using namespace std;
using namespace sts::likelihood;
using namespace sts::moves;
using namespace sts::particle;
using namespace sts::util;

const bpp::DNA DNA;
const bpp::RNA RNA;
const bpp::ProteicAlphabet AA;

std::vector<std::string> get_model_names()
{
    std::vector<string> models;
    // Nucleotide
    models.push_back("JCnuc");
    models.push_back("HKY");
    models.push_back("GTR");
    models.push_back("TN93");

    // Protein
    models.push_back("JTT");
    models.push_back("WAG");
    return models;
}

std::vector<std::string> get_bl_prop_dens()
{
    std::vector<string> prop_dens;
    prop_dens.push_back("expon");
    prop_dens.push_back("gamma");
    prop_dens.push_back("delta");
    prop_dens.push_back("unif2");
    return prop_dens;
}

/// Get the alphabet & substitution model associated with a name.

/// Model should match option from model_name_arg
shared_ptr<bpp::SubstitutionModel> model_for_name(const string name)
{
    shared_ptr<bpp::SubstitutionModel> model;
    // Nucleotide models
    if(name == "JCnuc" || name == "HKY" || name == "GTR" || name == "TN93") {
        if(name == "JCnuc") model = make_shared<bpp::JCnuc>(&DNA);
        else if(name == "HKY") model = make_shared<bpp::HKY85>(&DNA);
        else if(name == "GTR") model = make_shared<bpp::GTR>(&DNA);
        else if(name == "TN93") model = make_shared<bpp::TN93>(&DNA);
        else assert(false);
    } else {
        // Protein model
        if(name == "JTT") model = make_shared<bpp::JTT92>(&AA);
        else if(name == "WAG") model = make_shared<bpp::WAG01>(&AA);
        else assert(false);
    }
    return model;
}

int main(int argc, char** argv)
{
    TCLAP::CmdLine cmd("runs sts", ' ', STRINGIFY(STS_VERSION));

    TCLAP::UnlabeledValueArg<string> alignment(
        "alignment", "Input fasta alignment", true, "", "fasta alignment", cmd);
    TCLAP::ValueArg<string> output_path(
        "o", "out", "Where to write the output trees", false, "-", "tsv file", cmd);
    vector<string> all_models = get_model_names();
    TCLAP::ValuesConstraint<string> allowed_models(all_models);
    TCLAP::ValueArg<string> model_name(
        "m", "model-name", "Which substitution model to use", false, "JCnuc", &allowed_models, cmd);
    TCLAP::ValueArg<long> particle_count(
        "p", "particle-count", "Number of particles in the SMC", false, 1000, "#", cmd);
    TCLAP::SwitchArg no_compress("", "no-compress", "Do not compress the alignment to unique sites", cmd, false);
    TCLAP::ValueArg<int> bl_opt_steps(
        "", "bl-opt-steps", "Number of branch length optimization steps", false, 0, "#", cmd);
    std::vector<std::string> all_bl_prop_dens = get_bl_prop_dens();
    TCLAP::ValuesConstraint<string> allowed_bl_prop_dens(all_bl_prop_dens);
    TCLAP::ValueArg<string> bl_prop_dens(
        "", "bl-prop-dens", "Branch length proposal density", false, "expon", &allowed_bl_prop_dens, cmd);

    try {
        cmd.parse(argc, argv);
    } catch(TCLAP::ArgException &e) {
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
        return 1;
    }

    const long population_size = particle_count.getValue();
    ifstream in(alignment.getValue().c_str());
    string output_filename = output_path.getValue();
    ostream *output_stream;
    ofstream output_ofstream;
    if(output_filename == "-") {
        output_stream = &cout;
    } else {
        output_ofstream.open(output_filename);
        output_stream = &output_ofstream;
    }

    shared_ptr<bpp::SubstitutionModel> model = model_for_name(model_name.getValue());
    shared_ptr<bpp::SiteContainer> input_alignment = shared_ptr<bpp::SiteContainer>(read_alignment(in, model->getAlphabet()));
    shared_ptr<bpp::SiteContainer> aln = input_alignment;

    // Compress sites
    if(!no_compress.getValue())
        aln.reset(unique_sites(*aln));
    const int num_iters = aln->getNumberOfSequences();

    // Leaves
    vector<node> leaf_nodes;

    shared_ptr<online_calculator> calc = make_shared<online_calculator>();
    calc->initialize(aln, model);
    if(!no_compress.getValue())
        calc->set_weights(compressed_site_weights(*input_alignment, *aln));
    leaf_nodes.resize(num_iters);
    unordered_map<node, string> name_map;
    for(int i = 0; i < num_iters; i++) {
        leaf_nodes[i] = make_shared<phylo_node>(calc);
        calc->register_leaf(leaf_nodes[i], aln->getSequencesNames()[i]);
        name_map[leaf_nodes[i]] = aln->getSequencesNames()[i];
    }
    forest_likelihood fl(calc, leaf_nodes);

    rooted_merge::bl_proposal_fn chosen_bl_proposer, chosen_eb_bl_proposer;
    string bl_prop_dens_str = bl_prop_dens.getValue();
    if(bl_prop_dens_str == "expon") { // The exponential distribution with the supplied mean.
        auto loc_blp = exponential_branch_length_proposer(1.0);
        chosen_bl_proposer = loc_blp;
        chosen_eb_bl_proposer =
            eb_bl_proposer<exponential_branch_length_proposer>(fl, loc_blp, bl_opt_steps.getValue());
    }
    else if(bl_prop_dens_str == "gamma") { // The gamma distribution with shape = 2 with the supplied mean.
        auto loc_blp = gamma_branch_length_proposer(1.0);
        chosen_bl_proposer = loc_blp;
        chosen_eb_bl_proposer =
            eb_bl_proposer<gamma_branch_length_proposer>(fl, loc_blp, bl_opt_steps.getValue());
    }
    else if(bl_prop_dens_str == "delta") { // The delta distribution at the given mean.
        auto loc_blp = delta_branch_length_proposer(1.0);
        chosen_bl_proposer = loc_blp;
        chosen_eb_bl_proposer =
            eb_bl_proposer<delta_branch_length_proposer>(fl, loc_blp, bl_opt_steps.getValue());
    }
    else if(bl_prop_dens_str == "unif2") { // The uniform distribution on [0,2].
        auto loc_blp = uniform_branch_length_proposer(1.0); // The mean of the uniform distribution on [0,2] is 1.
        chosen_bl_proposer = loc_blp;
        chosen_eb_bl_proposer =
            eb_bl_proposer<uniform_branch_length_proposer>(fl, loc_blp, bl_opt_steps.getValue());
    }
    else {
        assert(false);
    }
    rooted_merge::bl_proposal_fn blp;
    if(!bl_opt_steps.getValue()){
        blp = chosen_bl_proposer;
    } else {
        blp = chosen_eb_bl_proposer;
    }

    rooted_merge smc_mv(fl, blp);
    smc_init init(fl);
    uniform_bl_mcmc_move mcmc_mv(fl, 0.1);

    try {

        // Initialise and run the sampler
        smc::sampler<particle> Sampler(population_size, SMC_HISTORY_NONE);
        smc::moveset<particle> Moveset(init, smc_mv, mcmc_mv);

        Sampler.SetResampleParams(SMC_RESAMPLE_STRATIFIED, 0.99);
        Sampler.SetMoveSet(Moveset);
        Sampler.Initialise();

        for(int n = 1 ; n < num_iters ; ++n) {
            const double ess = Sampler.IterateEss();

            double max_ll = -std::numeric_limits<double>::max();
            for(int i = 0; i < population_size; i++) {
                particle X = Sampler.GetParticleValue(i);
                // write the log likelihood
                double ll = fl(X);
                max_ll = max_ll > ll ? max_ll : ll;
            }
            cerr << "Iter " << n << " max ll " << max_ll << " ESS: " << ess << endl;
        }

        for(int i = 0; i < population_size; i++) {
            particle X = Sampler.GetParticleValue(i);
            // write the log likelihood
            *output_stream << fl(X) << "\t";
            // write out the tree under this particle
            write_tree(*output_stream, X->node, name_map);
        }
    }

    catch(smc::exception  e) {
        cerr << e;
        return e.lCode;
    }
    return 0;
}
