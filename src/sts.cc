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

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/shared_ptr.hpp>

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
boost::shared_ptr<bpp::SubstitutionModel> model_for_name(const string name)
{
    boost::shared_ptr<bpp::SubstitutionModel> model;
    // Nucleotide models
    if(name == "JCnuc" || name == "HKY" || name == "GTR" || name == "TN93") {
        if(name == "JCnuc") model = boost::make_shared<bpp::JCnuc>(&DNA);
        else if(name == "HKY") model = boost::make_shared<bpp::HKY85>(&DNA);
        else if(name == "GTR") model = boost::make_shared<bpp::GTR>(&DNA);
        else if(name == "TN93") model = boost::make_shared<bpp::TN93>(&DNA);
        else assert(false);
    } else {
        // Protein model
        if(name == "JTT") model = boost::make_shared<bpp::JTT92>(&AA);
        else if(name == "WAG") model = boost::make_shared<bpp::WAG01>(&AA);
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

template<class T>
size_t hash_value(const boost::shared_ptr<T>& aSptr)
{
};

/*
template<typename T>
struct MyClassEqual
{
  inline bool operator()(const MyClassPtr & p, const MyClassPtr & q)
  {
    return false; // implement
  }
};
*/
namespace boost {
namespace serialization {

template<class Archive>
void serialize(Archive & ar, smc::sampler<particle>& sampler, const unsigned int version)
{
    // oh man this is ugly.
    // but the particle array is private, so doing this temporarily to avoid changing smctc until we
    // decide on whether boost is useful.
    vector<particle> parts;
    for(int i=0; i<sampler.GetNumber(); i++){
        parts.push_back(sampler.GetParticleValue(i));
    }
    ar & boost::serialization::make_nvp("particles", parts);
}

template<class Archive>
void serialize(Archive & ar, sts::particle::phylo_particle& p, const unsigned int version)
{
    ar & boost::serialization::make_nvp("node",p.node);
    ar & boost::serialization::make_nvp("predecessor",p.predecessor);
}

template<class Archive>
void serialize(Archive & ar, sts::particle::phylo_node& n, const unsigned int version)
{
    ar & boost::serialization::make_nvp("height",n.height);
    ar & boost::serialization::make_nvp("child1",n.child1);
    ar & boost::serialization::make_nvp("child2",n.child2);
}

template<class Archive>
void serialize(Archive & ar, sts::particle::edge& e, const unsigned int version)
{
    ar & boost::serialization::make_nvp("length",e.length);
    ar & boost::serialization::make_nvp("node",e.node);
}

} // namespace serialization
} // namespace boost


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
    if (output_filename == "-") {
        output_stream = &cout;
    } else {
        output_ofstream.open(output_filename);
        output_stream = &output_ofstream;
    }

    boost::shared_ptr<bpp::SubstitutionModel> model = model_for_name(model_name.getValue());
    boost::shared_ptr<bpp::SiteContainer> aln = boost::shared_ptr<bpp::SiteContainer>(read_alignment(in, model->getAlphabet()));

    // Compress sites
    if(!no_compress.getValue())
        aln.reset(unique_sites(*aln));
    const int num_iters = aln->getNumberOfSequences();

    // Leaves
    vector<node> leaf_nodes;

    boost::shared_ptr<online_calculator> calc = boost::make_shared<online_calculator>();
    calc->initialize(aln, model);
    leaf_nodes.resize(num_iters);
    unordered_map<node, string> name_map;
    for(int i = 0; i < num_iters; i++) {
        leaf_nodes[i] = make_shared<phylo_node>(calc);
        calc->register_node(leaf_nodes[i]);
        name_map[leaf_nodes[i]] = aln->getSequencesNames()[i];
    }
    forest_likelihood fl(calc, leaf_nodes);
    rooted_merge smc_mv(fl);
    smc_init init(fl);
    uniform_bl_mcmc_move mcmc_mv(fl, 0.1);
    
    ofstream boost_out("boost.out");

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
            }
            cerr << "Iter " << n << " max ll " << max_ll << endl;
//            boost::archive::text_oarchive oa(boost_out);
//            oa << Sampler;
            boost::archive::xml_oarchive oa(boost_out);
            oa << boost::serialization::make_nvp("sampler",Sampler);
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
