#include "log/sampler.h"
#include "log/json_logger.h"
#include "forest_likelihood.h"
#include "online_calculator.h"
#include "node.h"
#include "state.h"
#include "util.h"
#include "eb_bl_proposer.h"
#include "uniform_bl_mcmc_move.h"
#include "rooted_merge.h"
#include "smc_init.h"

#include "delta_branch_length_proposer.h"
#include "exponential_branch_length_proposer.h"
#include "gamma_branch_length_proposer.h"
#include "uniform_branch_length_proposer.h"


#include <Bpp/Phyl/Model/GTR.h>
#include <Bpp/Phyl/Model/HKY85.h>
#include <Bpp/Phyl/Model/JCnuc.h>
#include <Bpp/Phyl/Model/JTT92.h>
#include <Bpp/Phyl/Model/TN93.h>
#include <Bpp/Phyl/Model/WAG01.h>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Alphabet/RNA.h>
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>
#include <Bpp/Seq/Container/SiteContainer.h>

#include "tclap/CmdLine.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <sstream>
#include <stack>
#include <unordered_map>

#define _STRINGIFY(s) #s
#define STRINGIFY(s) _STRINGIFY(s)

using namespace std;
using namespace sts::log;
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

std::vector<std::string> get_bl_dens()
{
    std::vector<string> dens;
    dens.push_back("expon");
    dens.push_back("gamma");
    dens.push_back("delta");
    dens.push_back("unif2");
    return dens;
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
    TCLAP::ValueArg<string> log_path(
        "l", "log", "Filename for particle state log (JSON format)", false, "", "json file", cmd);
    vector<string> all_models = get_model_names();
    TCLAP::ValuesConstraint<string> allowed_models(all_models);
    TCLAP::ValueArg<string> model_name(
        "m", "model-name", "Which substitution model to use", false, "JCnuc", &allowed_models, cmd);
    TCLAP::ValueArg<long> particle_count(
        "p", "particle-count", "Number of particles in the SMC", false, 1000, "#", cmd);
    TCLAP::SwitchArg no_compress("", "no-compress", "Do not compress the alignment to unique sites", cmd, false);
    TCLAP::ValueArg<int> bl_opt_steps(
        "", "bl-opt-steps", "Number of branch length optimization steps", false, 0, "#", cmd);
    std::vector<std::string> all_bl_dens = get_bl_dens();
    TCLAP::ValuesConstraint<string> allowed_bl_dens(all_bl_dens);
    TCLAP::ValueArg<string> bl_dens(
        "", "bl-dens", "Branch length prior & proposal density", false, "expon", &allowed_bl_dens, cmd);
    TCLAP::SwitchArg verify_ll("", "verify-cached-ll", "Verify cached log-likelihoods", cmd, false);

    try {
        cmd.parse(argc, argv);
    } catch(TCLAP::ArgException &e) {
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
        return 1;
    }

    const long population_size = particle_count.getValue();
    ifstream in(alignment.getValue());
    if(in.bad()) {
        cerr << "Cannot read from " << alignment.getValue() << endl;
        return 1;
    }
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
    vector<Node_ptr> leaf_nodes;

    shared_ptr<Online_calculator> calc = make_shared<Online_calculator>();
    calc->verify_cached_ll = verify_ll.getValue();
    calc->initialize(aln, model);
    if(!no_compress.getValue())
        calc->set_weights(compressed_site_weights(*input_alignment, *aln));
    leaf_nodes.resize(num_iters);
    unordered_map<Node_ptr, string> node_name_map;
    for(int i = 0; i < num_iters; i++) {
        leaf_nodes[i] = make_shared<Node>(calc);
        calc->register_leaf(leaf_nodes[i], aln->getSequencesNames()[i]);
        node_name_map[leaf_nodes[i]] = aln->getSequencesNames()[i];
    }
    Forest_likelihood forest_likelihood(calc, leaf_nodes);

    Rooted_merge::Bl_proposal_fn chosen_bl_proposer, chosen_eb_bl_proposer;
    string bl_dens_str = bl_dens.getValue();
    if(bl_dens_str == "expon") { // The exponential distribution with the supplied mean.
        auto loc_blp = Exponential_branch_length_proposer(1.0);
        chosen_bl_proposer = loc_blp;
        chosen_eb_bl_proposer =
            Eb_bl_proposer<Exponential_branch_length_proposer>(forest_likelihood, loc_blp, bl_opt_steps.getValue());
    } else if(bl_dens_str == "gamma") { // The gamma distribution with shape = 2 with the supplied mean.
        auto loc_blp = Gamma_branch_length_proposer(1.0);
        chosen_bl_proposer = loc_blp;
        chosen_eb_bl_proposer =
            Eb_bl_proposer<Gamma_branch_length_proposer>(forest_likelihood, loc_blp, bl_opt_steps.getValue());
    } else if(bl_dens_str == "delta") { // The delta distribution at the given mean.
        auto loc_blp = Delta_branch_length_proposer(1.0);
        chosen_bl_proposer = loc_blp;
        chosen_eb_bl_proposer =
            Eb_bl_proposer<Delta_branch_length_proposer>(forest_likelihood, loc_blp, bl_opt_steps.getValue());
    } else if(bl_dens_str == "unif2") { // The uniform distribution on [0,2].
        auto loc_blp = Uniform_branch_length_proposer(1.0); // The mean of the uniform distribution on [0,2] is 1.
        chosen_bl_proposer = loc_blp;
        chosen_eb_bl_proposer =
            Eb_bl_proposer<Uniform_branch_length_proposer>(forest_likelihood, loc_blp, bl_opt_steps.getValue());
    } else {
        assert(false);
    }
    Rooted_merge::Bl_proposal_fn final_bl_proposer;
    if(!bl_opt_steps.getValue()) {
        final_bl_proposer = chosen_bl_proposer;
    } else {
        final_bl_proposer = chosen_eb_bl_proposer;
    }

    Rooted_merge smc_mv(forest_likelihood, final_bl_proposer);
    Smc_init init(forest_likelihood);
    Uniform_bl_mcmc_move mcmc_mv(forest_likelihood, 0.1);

    ofstream json_out;
    unique_ptr<Json_logger> logger;
    if(!log_path.getValue().empty()) {
        json_out.open(log_path.getValue());
        logger.reset(new Json_logger(json_out));
    }

    try {

        // Initialize and run the sampler.
        smc::sampler<Particle> Sampler(population_size, SMC_HISTORY_NONE);
        smc::moveset<Particle> Moveset(init, smc_mv, mcmc_mv);

        Sampler.SetResampleParams(SMC_RESAMPLE_STRATIFIED, 0.99);
        Sampler.SetMoveSet(Moveset);
        Sampler.Initialise();

        for(int n = 1 ; n < num_iters ; ++n) {
            const double ess = Sampler.IterateEss();

            double max_ll = -std::numeric_limits<double>::max();
            for(int i = 0; i < population_size; i++) {
                Particle X = Sampler.GetParticleValue(i);
                double ll = forest_likelihood(X);
                max_ll = max_ll > ll ? max_ll : ll;
            }
            cerr << "Iter " << n << " max ll " << max_ll << " ESS: " << ess << endl;
            if(logger) logger->log(Sampler, node_name_map);
        }

        if(logger) logger->write();

        for(int i = 0; i < population_size; i++) {
            Particle X = Sampler.GetParticleValue(i);
            // Write the log likelihood.
            *output_stream << forest_likelihood(X) << "\t";
            // Write out the tree under this particle.
            write_tree(*output_stream, X->node, node_name_map);
        }
    }

    catch(smc::exception  e) {
        cerr << e;
        return e.lCode;
    }
    return 0;
}
