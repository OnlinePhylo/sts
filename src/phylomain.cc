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
#include <Bpp/Phyl/Model/JCnuc.h>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Io/IoSequenceFactory.h>
#include <Bpp/Seq/Io/ISequence.h>

using namespace std;

bpp::SiteContainer* read_alignment(istream &in, bpp::Alphabet *alphabet)
{
    bpp::Fasta r;
    bpp::SiteContainer *sequences = new bpp::VectorSiteContainer(alphabet);
    bpp::Sequence *seq = new bpp::BasicSequence(alphabet);

    while(r.nextSequence(in, *seq)) {
        sequences->addSequence(*seq, true);
        seq = new bpp::BasicSequence(alphabet);
    }

    // One more seq allocated
    delete seq;

    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sequences);
    cerr << sequences->getNumberOfSequences() << " sequences" << endl << sequences->getNumberOfSites() << " sites" << endl;

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
        if(!visited_id(visited, cur->child1->id)) {
            out << "(";
            s.push(cur->child1);
            continue;
        } else if(!visited_id(visited, cur->child2->id)) {
            out << ":" << cur->dist1 << ",";
            s.push(cur->child2);
            continue;
        }
        out << ":" << cur->dist2 << ")";
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
            if(node_x.find(cur->child1->id) == node_x.end()) {
                s.push(cur->child1);
                continue;
            } else if(node_x.find(cur->child2->id) == node_x.end()) {
                s.push(cur->child2);
                continue;
            }
            node_x[cur->id] = (node_x[cur->child1->id] + node_x[cur->child2->id]) / 2.0;
            node_y[cur->id] = margin + cur->height * height_scaler;
            left.push_back(node_x[cur->child1->id]);
            right.push_back(node_x[cur->child2->id]);
            lfrom.push_back(node_y[cur->child1->id]);
            rfrom.push_back(node_y[cur->child2->id]);
            to.push_back(node_y[cur->id]);
            s.pop();
        }
    }
    // Write out a list of the particle drawing instructions
    for(int i = 0; i < left.size(); i++) {
        out << left[i] << "\t" << right[i] << "\t" << lfrom[i] << "\t" << rfrom[i] << "\t" << to[i] << endl;
    }
}

int main(int argc, char** argv)
{
    if(argc != 2) {
        cerr << "Usage: phylo <fasta alignment>\n\n";
        return -1;
    }
    long population_size = 1000;

    string file_name = argv[1];

    ifstream in(file_name.c_str());
    bpp::DNA dna;
    model.reset(new bpp::JCnuc(&dna));
    aln.reset(read_alignment(in, &dna));

    ofstream viz_pipe("viz_data.csv");

    long lIterates = aln->getNumberOfSequences();

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


