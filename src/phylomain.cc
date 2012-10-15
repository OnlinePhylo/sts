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

using namespace std;

void read_alignment(istream& in, vector< pair< string, string > >& aln)
{
    string line, cur_seq, name;
    while(getline(in, line)) {
        if(line[0] == '>') {
            if(cur_seq.size() > 0) {
                aln.push_back(make_pair(name, cur_seq));
            }
            name = line.substr(1);
            cur_seq = "";
        } else if(line[0] == '#') {
            break;
        } else {
            cur_seq += line;
        }
    }
    aln.push_back(make_pair(name, cur_seq));
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
    vector< bool > visited;
    stack< shared_ptr< phylo_node > > s;
    s.push(root);
    while(s.size() > 0) {
        shared_ptr< phylo_node > cur = s.top();
        if(cur->child1 == NULL) {
            out << aln[leaf_sequence_ids[cur]].first;
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
    int leaf_unit = (viz_width - 2 * margin) / aln.size();
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

template <typename T>
class phylo_sampler : public smc::sampler< T >
{
public:
    phylo_sampler(int population_size, int rng_seed)
     : smc::sampler< T >( population_size, SMC_HISTORY_NONE, gsl_rng_default, rng_seed ) 
    {
        // Initialise the sampler
        smc::moveset< T > Moveset(fInitialise, fMove, fMoveNodeAgeMCMC);
        this->SetResampleParams(SMC_RESAMPLE_STRATIFIED, 0.99);
        this->SetMoveSet(Moveset);
        this->Initialise();

        viz_pipe.open("viz_data.csv");
        seq_count = 0;
    }

    void process_new_data()
    {
        std::vector< std::string > just_the_seqs_maam;
        for(int i = seq_count; i < aln.size(); i++)
            just_the_seqs_maam.push_back(aln[i].second);
        std::vector< shared_ptr< phylo_node > > new_leaves;
        calc.initialize(just_the_seqs_maam, new_leaves);
        leaf_nodes.insert( leaf_nodes.end(), new_leaves.begin(), new_leaves.end() );        
        for(int i = 0; i < new_leaves.size(); i++) {
            leaf_sequence_ids[new_leaves[i]]=seq_count+i;
        }

    try {
        
        int lIterates = aln.size() - seq_count;
        cerr << "Adding " << lIterates << " new sequences\n";
        for(int n = 1 ; n < lIterates ; ++n) {
            this->Iterate();

            double max_ll = -std::numeric_limits<double>::max();
            for(int i = 0; i < this->GetNumber(); i++) {
                particle X = this->GetParticleValue(i);
                // write the log likelihood
                double ll = logLikelihood(lIterates, X);
                max_ll = max_ll > ll ? max_ll : ll;

                write_forest_viz(viz_pipe, X.pp);
            }
            viz_pipe << "############## End of generation ##############\n";
            cerr << "Iter " << seq_count + n << " max ll " << max_ll << endl;
        }

        for(int i = 0; i < this->GetNumber(); i++) {
            particle X = this->GetParticleValue(i);
            // write the log likelihood
            cout << logLikelihood(lIterates, X) << "\t";
            // write out the tree under this particle
            write_tree(cout, X.pp->node);
        }

        seq_count = aln.size();
    }

    catch(smc::exception  e) {
        cerr << e;
        exit(e.lCode);
    }

    }
private:
    ofstream viz_pipe;
    int seq_count;
};


int main(int argc, char** argv)
{
    if(argc < 2) {
        cerr << "Usage: phylo <fasta alignment> [seed]\n\n";
        return -1;
    }
    long population_size = 3000;
    int seed = 666;
    if(argc == 3) seed = atoi(argv[2]);
    string file_name = argv[1];
    ifstream in(file_name.c_str());
    phylo_sampler< particle > sampler( population_size, seed );

    int cur_aln_size = 0;
    while(in.good()){
        read_alignment(in, aln);
        sampler.process_new_data();
    }
}


