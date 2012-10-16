#include "smctc.hh"
#include "phylofunc.hh"
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

bpp::SiteContainer* read_alignment(istream &in, bpp::Alphabet *alphabet)
{
    // Holy boilerplate - Bio++ won't allow reading FASTA files as alignments
    bpp::IOSequenceFactory fac;
    bpp::ISequence *reader = fac.createReader(bpp::IOSequenceFactory::FASTA_FORMAT);
    bpp::SequenceContainer *seqs = reader->read(in, alphabet);

    // Have to look up by name
    vector<string> names = seqs->getSequencesNames();
    bpp::SiteContainer *sequences = new bpp::VectorSiteContainer(alphabet);

    for(auto it = names.begin(), end = names.end(); it != end; ++it) {
        sequences->addSequence(seqs->getSequence(*it), true);
    }

    // One more seq allocated
    delete seqs;

    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sequences);

    return sequences;
}

bool check_visited(vector< bool >& visited, int id)
{
    // ensure visited has enough space allocated to store the id
    // if not, resize it large enough and leave some wiggle to prevent frequent resizings
    if(id >= visited.size()) {
        visited.resize(id + 100);
    }
}

bool visited_id(vector< bool >& visited, int id)
{
    check_visited(visited, id);
    return visited[id];
}

bool set_visited_id(vector< bool >& visited, int id)
{
    check_visited(visited, id);
    visited[id] = true;
}

void write_tree(ostream& out, shared_ptr< phylo_node > root)
{
    vector<string> names = aln->getSequencesNames();
    vector< bool > visited;
    stack< shared_ptr< phylo_node > > s;
    s.push(root);
    while(s.size() > 0) {
        shared_ptr< phylo_node > cur = s.top();
        if(cur->child1 == NULL) {
            out << names[cur->id];
            set_visited_id(visited, cur->id);
            s.pop();
            continue;
        }
        if(!visited_id(visited, cur->child1->node->id)) {
            out << "(";
            s.push(cur->child1->node);
            continue;
        } else if(!visited_id(visited, cur->child2->node->id)) {
            out << ":" << cur->child1->length << ",";
            s.push(cur->child2->node);
            continue;
        }
        out << ":" << cur->child2->length << ")";
        set_visited_id(visited, cur->id);
        s.pop();
    }
    out << ";\n";
}

void write_forest_viz(ostream& out, shared_ptr< phylo_particle > part)
{
    int viz_width = 640;
    int viz_height = 320;
    int margin = 20;
    int leaf_unit = (viz_width - 2 * margin) / aln->getNumberOfSequences();
    float root_height_limit = 3.0;
    float height_scaler = (viz_height - margin * 2) / root_height_limit;
    vector< float > left;
    vector< float > right;
    vector< float > lfrom;
    vector< float > rfrom;
    vector< float > to;
    unordered_map< int, float > node_y;
    unordered_map< int, float > node_x;
    for(shared_ptr< phylo_particle > p = part; p != NULL; p = p->predecessor) {
        shared_ptr< phylo_node > root = p->node;
        if(root == NULL || node_x.find(root->id) != node_x.end())  continue;

        stack< shared_ptr< phylo_node > > s;
        s.push(root);
        while(s.size() > 0) {
            shared_ptr< phylo_node > cur = s.top();
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

std::vector<std::string> get_model_names() {
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
std::pair<bpp::Alphabet*, bpp::SubstitutionModel*> model_for_name(std::string name)
{
    bpp::SubstitutionModel *model;
    bpp::Alphabet *alphabet;
    // Nucleotide models
    if(name == "JCnuc" || name == "HKY" || name == "GTR" || name == "TN93") {
        bpp::DNA *dna = new bpp::DNA();
        alphabet = dna;
        if(name == "JCnuc") model = new bpp::JCnuc(dna);
        else if (name == "HKY") model = new bpp::HKY85(dna);
        else if (name == "GTR") model = new bpp::GTR(dna);
        else if (name == "TN93") model = new bpp::TN93(dna);
        else assert(false);
    } else {
        bpp::ProteicAlphabet *aa = new bpp::ProteicAlphabet();
        alphabet = aa;
        // Protein model
        if(name == "JTT") model = new bpp::JTT92(aa);
        else if(name == "WAG") model = new bpp::WAG01(aa);
        else assert(false);
    }
    return std::pair<bpp::Alphabet*,bpp::SubstitutionModel*>(alphabet, model);
}

int main(int argc, char** argv)
{
    TCLAP::CmdLine cmd("runs sts", ' ', STRINGIFY(STS_VERSION));

    TCLAP::UnlabeledValueArg<string> alignment(
        "alignment", "Input fasta alignment", true, "", "fasta alignment", cmd);
    vector<string> all_models = get_model_names();
    TCLAP::ValuesConstraint<string> allowed_models(all_models);
    TCLAP::ValueArg<string> model_name(
        "m", "model-name", "Which substitution model to use", false, "JC69", &allowed_models, cmd);
    TCLAP::ValueArg<long> particle_count(
            "p", "particle-count", "Number of particles in the SMC", false, 1000, "#", cmd);

    try {
        cmd.parse(argc, argv);
    } catch (TCLAP::ArgException &e) {
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
        return 1;
    }

    const long population_size = particle_count.getValue();
    ifstream in(alignment.getValue().c_str());


    string model_name_string = model_name.getValue();
    auto alpha_model = model_for_name(model_name_string);
    model.reset(alpha_model.second);

    aln.reset(read_alignment(in, alpha_model.first));

    ofstream viz_pipe("viz_data.csv");

    const long lIterates = aln->getNumberOfSequences();

    try {
        leaf_nodes.resize(lIterates);
        for(int i = 0; i < lIterates; i++) {
            leaf_nodes[i] = make_shared< phylo_node >();
            leaf_nodes[i]->id = i;
        }

        // Initialise and run the sampler
        smc::sampler<particle> Sampler(population_size, SMC_HISTORY_NONE);
        smc::moveset<particle> Moveset(fInitialise, fMove, fMoveNodeAgeMCMC);

        Sampler.SetResampleParams(SMC_RESAMPLE_STRATIFIED, 0.99);
        Sampler.SetMoveSet(Moveset);
        Sampler.Initialise();

        for(int n = 1 ; n < lIterates ; ++n) {
            Sampler.Iterate();

            double max_ll = -std::numeric_limits<double>::max();
            for(int i = 0; i < population_size; i++) {
                particle X = Sampler.GetParticleValue(i);
                // write the log likelihood
                double ll = logLikelihood(lIterates, X);
                max_ll = max_ll > ll ? max_ll : ll;

                write_forest_viz(viz_pipe, X.pp);
            }
            viz_pipe << "############## End of generation ##############\n";
            cerr << "Iter " << n << " max ll " << max_ll << endl;
        }

        for(int i = 0; i < population_size; i++) {
            particle X = Sampler.GetParticleValue(i);
            // write the log likelihood
            cout << logLikelihood(lIterates, X) << "\t";
            // write out the tree under this particle
            write_tree(cout, X.pp->node);
        }
    }

    catch(smc::exception  e) {
        cerr << e;
        exit(e.lCode);
    }
}
