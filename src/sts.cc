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
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>
#include <Bpp/Seq/Container/SequenceContainer.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Io/IoSequenceFactory.h>
#include <Bpp/Seq/Io/ISequence.h>
#include "tclap/CmdLine.h"

#define _STRINGIFY(s) #s
#define STRINGIFY(s) _STRINGIFY(s)

using namespace std;
using namespace sts::likelihood;
using namespace sts::moves;
using namespace sts::particle;

const bpp::DNA DNA;
const bpp::RNA RNA;
const bpp::ProteicAlphabet AA;

bpp::SiteContainer* read_alignment(istream &in, const bpp::Alphabet *alphabet)
{
    // Holy boilerplate - Bio++ won't allow reading FASTA files as alignments
    bpp::IOSequenceFactory fac;
    unique_ptr<bpp::ISequence> reader = unique_ptr<bpp::ISequence>(
                                            fac.createReader(bpp::IOSequenceFactory::FASTA_FORMAT));
    unique_ptr<bpp::SequenceContainer> seqs = unique_ptr<bpp::SequenceContainer>(reader->read(in, alphabet));

    // Have to look up by name
    vector<string> names = seqs->getSequencesNames();
    bpp::SiteContainer *sequences = new bpp::VectorSiteContainer(alphabet);

    for(auto it = names.begin(), end = names.end(); it != end; ++it) {
        sequences->addSequence(seqs->getSequence(*it), true);
    }

    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sequences);

    return sequences;
}

void write_forest_viz(ostream& out, const shared_ptr<phylo_particle> part, const size_t sequence_count)
{
    int viz_width = 640;
    int viz_height = 320;
    int margin = 20;
    int leaf_unit = (viz_width - 2 * margin) / sequence_count;
    float root_height_limit = 3.0;
    float height_scaler = (viz_height - margin * 2) / root_height_limit;
    vector<float> left;
    vector<float> right;
    vector<float> lfrom;
    vector<float> rfrom;
    vector<float> to;
    unordered_map<int, float> node_y;
    unordered_map<int, float> node_x;
    for(shared_ptr<phylo_particle> p = part; p != NULL; p = p->predecessor) {
        shared_ptr<phylo_node> root = p->node;
        if(root == NULL || node_x.find(root->id) != node_x.end())  continue;

        stack<shared_ptr<phylo_node>> s;
        s.push(root);
        while(s.size() > 0) {
            shared_ptr<phylo_node> cur = s.top();
            if(cur->child1 == NULL) {
                node_y[cur->id] = margin;
                node_x[cur->id] = margin + leaf_unit * cur->id;
                s.pop();
                continue;
            }
            if(node_x.find(cur->child1->node->id) == node_x.end()) {
                s.push(cur->child1->node);
                continue;
            } else if(node_x.find(cur->child2->node->id) == node_x.end()) {
                s.push(cur->child2->node);
                continue;
            }
            node_x[cur->id] = (node_x[cur->child1->node->id] + node_x[cur->child2->node->id]) / 2.0;
            node_y[cur->id] = margin + cur->height * height_scaler;
            left.push_back(node_x[cur->child1->node->id]);
            right.push_back(node_x[cur->child2->node->id]);
            lfrom.push_back(node_y[cur->child1->node->id]);
            rfrom.push_back(node_y[cur->child2->node->id]);
            to.push_back(node_y[cur->id]);
            s.pop();
        }
    }
    // Write out a list of the particle drawing instructions
    for(int i = 0; i < left.size(); i++) {
        out << left[i] << "\t" << right[i] << "\t" << lfrom[i] << "\t" << rfrom[i] << "\t" << to[i] << endl;
    }
}

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

/// Get the alphabet & substitution model associated with a name.

/// Model should match option from model_name_arg
std::shared_ptr<bpp::SubstitutionModel> model_for_name(const std::string name)
{
    std::shared_ptr<bpp::SubstitutionModel> model;
    // Nucleotide models
    if(name == "JCnuc" || name == "HKY" || name == "GTR" || name == "TN93") {
        if(name == "JCnuc") model = std::make_shared<bpp::JCnuc>(&DNA);
        else if(name == "HKY") model = std::make_shared<bpp::HKY85>(&DNA);
        else if(name == "GTR") model = std::make_shared<bpp::GTR>(&DNA);
        else if(name == "TN93") model = std::make_shared<bpp::TN93>(&DNA);
        else assert(false);
    } else {
        // Protein model
        if(name == "JTT") model = std::make_shared<bpp::JTT92>(&AA);
        else if(name == "WAG") model = std::make_shared<bpp::WAG01>(&AA);
        else assert(false);
    }
    return model;
}

static bpp::SiteContainer* unique_sites(const bpp::SiteContainer& sites)
{
    bpp::SiteContainer *compressed = bpp::PatternTools::shrinkSiteSet(sites);

    if(compressed->getNumberOfSites() < sites.getNumberOfSites())
        cerr << "Reduced from "
             << sites.getNumberOfSites()
             << " to " << compressed->getNumberOfSites()
             << " sites"
             << endl;

    return compressed;
}

int main(int argc, char** argv)
{
    TCLAP::CmdLine cmd("runs sts", ' ', STRINGIFY(STS_VERSION));

    TCLAP::UnlabeledValueArg<string> alignment(
        "alignment", "Input fasta alignment", true, "", "fasta alignment", cmd);
    vector<string> all_models = get_model_names();
    TCLAP::ValuesConstraint<string> allowed_models(all_models);
    TCLAP::ValueArg<string> model_name(
        "m", "model-name", "Which substitution model to use", false, "JCnuc", &allowed_models, cmd);
    TCLAP::ValueArg<long> particle_count(
        "p", "particle-count", "Number of particles in the SMC", false, 1000, "#", cmd);
    TCLAP::SwitchArg no_compress("", "no-compress", "Do not compress the alignment to unique sites", cmd, false);

    try {
        cmd.parse(argc, argv);
    } catch(TCLAP::ArgException &e) {
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
        return 1;
    }

    const long population_size = particle_count.getValue();
    ifstream in(alignment.getValue().c_str());

    std::shared_ptr<bpp::SubstitutionModel> model = model_for_name(model_name.getValue());
    std::shared_ptr<bpp::SiteContainer> aln = std::shared_ptr<bpp::SiteContainer>(read_alignment(in, model->getAlphabet()));

    // Compress sites
    if(!no_compress.getValue())
        aln.reset(unique_sites(*aln));
    const int num_iters = aln->getNumberOfSequences();

    // Leaves
    std::vector<std::shared_ptr<phylo_node>> leaf_nodes;

    ofstream viz_pipe("viz_data.csv");

    std::shared_ptr<online_calculator> calc = std::make_shared<online_calculator>();
    leaf_nodes.resize(num_iters);
    for(int i = 0; i < num_iters; i++) {
        leaf_nodes[i] = make_shared<phylo_node>(calc);
        leaf_nodes[i]->id = i;
    }
    calc->initialize(aln, model);
    forest_likelihood fl(calc, leaf_nodes);

    // ML Optimization test
    exponential_branch_length_proposer exp_prop(1.0);
    eb_bl_proposer<exponential_branch_length_proposer> p(fl, exp_prop, 4);

    rooted_merge smc_mv(fl, p);
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
            Sampler.Iterate();

            double max_ll = -std::numeric_limits<double>::max();
            for(int i = 0; i < population_size; i++) {
                particle X = Sampler.GetParticleValue(i);
                // write the log likelihood
                double ll = fl(X);
                max_ll = max_ll > ll ? max_ll : ll;

                write_forest_viz(viz_pipe, X, aln->getNumberOfSequences());
            }
            viz_pipe << "############## End of generation ##############\n";
            cerr << "Iter " << n << " max ll " << max_ll << endl;
        }

        for(int i = 0; i < population_size; i++) {
            particle X = Sampler.GetParticleValue(i);
            // write the log likelihood
            cout << fl(X) << "\t";
            // write out the tree under this particle
            write_tree(cout, X->node, aln->getSequencesNames());
        }
    }

    catch(smc::exception  e) {
        cerr << e;
        return e.lCode;
    }
    return 0;
}
