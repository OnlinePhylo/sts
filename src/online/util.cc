/// \file util.cc
/// \brief Implementation of utility functions.
///
#include "util.h"

#ifndef NO_BEAGLE
#include "libhmsbeagle/beagle.h"
#endif

#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Seq/Container/SequenceContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Io/IoSequenceFactory.h>
#include <Bpp/Seq/Io/ISequence.h>
#include <Bpp/Seq/SiteTools.h>

#include <cassert>
#include <memory>
#include <numeric>
#include <stack>
#include <stdexcept>
#include <unordered_set>
#include <utility>

/// STS namespace
namespace sts
{
/// \brief Utility functions
/// \author Metatangle, inc.
namespace util
{


/// Read an alignment from a stream
/// \param in Input stream
/// \param alphabet The alphabet to use.
bpp::SiteContainer* read_alignment(std::istream &in, const bpp::Alphabet *alphabet)
{
    bpp::IoSequenceFactory fac;
    std::unique_ptr<bpp::ISequence> reader = std::unique_ptr<bpp::ISequence>(
                fac.createReader(bpp::IoSequenceFactory::FASTA_FORMAT));
    std::unique_ptr<bpp::SequenceContainer> raw_seqs(reader->readSequences(in, alphabet));
    std::unique_ptr<bpp::SiteContainer> sequences(new bpp::VectorSiteContainer(*raw_seqs));
    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sequences);

    return sequences.release();
}

/// Determine new site weights after compression.
/// \param orig Original sites
/// \param compressed sites after compression to unique sites
/// \returns A vector, where the value at each position is the appropriate weight for <c>compressed</c> site \c i
std::vector<double> compressed_site_weights(const bpp::SiteContainer& orig, const bpp::SiteContainer& compressed)
{
    assert(compressed.getNumberOfSites() <= orig.getNumberOfSites());

    std::vector<double> result(compressed.getNumberOfSites(), 0);

    // Get the first index at which each site from orig appears in compressed.
    std::vector<int> m = bpp::PatternTools::getIndexes(orig, compressed);
    for(unsigned int i = 0; i < m.size(); ++i) {
        assert(m[i] >= 0);
        ++result[m[i]];
    }

    unsigned int tot = std::accumulate(result.begin(), result.end(), 0);
    assert(tot == orig.getNumberOfSites());

    return result;
}

/// Get the unique sites in an alignment
///
/// \param sites Original sites
/// \param verbose Print a message if sites are compressed?
/// \returns A sequence container containing only the unique sites from \c sites.
bpp::SiteContainer* unique_sites(const bpp::SiteContainer& sites, bool verbose)
{
    bpp::SiteContainer *compressed = bpp::PatternTools::shrinkSiteSet(sites);

    if(verbose && compressed->getNumberOfSites() < sites.getNumberOfSites())
        std::clog << "Reduced from "
                  << sites.getNumberOfSites()
                  << " to " << compressed->getNumberOfSites()
                  << " sites"
                  << std::endl;

    return compressed;
}

#ifndef NO_BEAGLE
std::string beagle_errstring(const int beagle_error_code)
{
    switch(beagle_error_code) {
        case BEAGLE_SUCCESS:                      return "BEAGLE_SUCCESS";
        case BEAGLE_ERROR_GENERAL:                return "BEAGLE_ERROR_GENERAL";
        case BEAGLE_ERROR_OUT_OF_MEMORY:          return "BEAGLE_ERROR_OUT_OF_MEMORY";
        case BEAGLE_ERROR_UNIDENTIFIED_EXCEPTION: return "BEAGLE_ERROR_UNIDENTIFIED_EXCEPTION";
        case BEAGLE_ERROR_UNINITIALIZED_INSTANCE: return "BEAGLE_ERROR_UNINITIALIZED_INSTANCE";
        case BEAGLE_ERROR_OUT_OF_RANGE:           return "BEAGLE_ERROR_OUT_OF_RANGE";
        case BEAGLE_ERROR_NO_RESOURCE:            return "BEAGLE_ERROR_NO_RESOURCE";
        case BEAGLE_ERROR_NO_IMPLEMENTATION:      return "BEAGLE_ERROR_NO_IMPLEMENTATION";
        case BEAGLE_ERROR_FLOATING_POINT:         return "BEAGLE_ERROR_FLOATING_POINT";
        default: return "Unknown Beagle error: " + std::to_string(beagle_error_code);
    }
}

void beagle_check(int return_code)
{
    if(return_code != BEAGLE_SUCCESS)
        throw std::runtime_error(sts::util::beagle_errstring(return_code));
}

// Fucntion from MrBayes
void beagle_print_flags(long inFlags)
{
	char *names[] = { "PROCESSOR_CPU",
		"PROCESSOR_GPU",
		"PROCESSOR_FPGA",
		"PROCESSOR_CELL",
		"PRECISION_DOUBLE",
		"PRECISION_SINGLE",
		"COMPUTATION_ASYNCH",
		"COMPUTATION_SYNCH",
		"EIGEN_REAL",
		"EIGEN_COMPLEX",
		"SCALING_MANUAL",
		"SCALING_AUTO",
		"SCALING_ALWAYS",
		"SCALING_DYNAMIC",
		"SCALERS_RAW",
		"SCALERS_LOG",
		"VECTOR_NONE",
		"VECTOR_SSE",
		"THREADING_NONE",
		"THREADING_OPENMP"
	};
	long flags[] = { BEAGLE_FLAG_PROCESSOR_CPU,
		BEAGLE_FLAG_PROCESSOR_GPU,
		BEAGLE_FLAG_PROCESSOR_FPGA,
		BEAGLE_FLAG_PROCESSOR_CELL,
		BEAGLE_FLAG_PRECISION_DOUBLE,
		BEAGLE_FLAG_PRECISION_SINGLE,
		BEAGLE_FLAG_COMPUTATION_ASYNCH,
		BEAGLE_FLAG_COMPUTATION_SYNCH,
		BEAGLE_FLAG_EIGEN_REAL,
		BEAGLE_FLAG_EIGEN_COMPLEX,
		BEAGLE_FLAG_SCALING_MANUAL,
		BEAGLE_FLAG_SCALING_AUTO,
		BEAGLE_FLAG_SCALING_ALWAYS,
		BEAGLE_FLAG_SCALING_DYNAMIC,
		BEAGLE_FLAG_SCALERS_RAW,
		BEAGLE_FLAG_SCALERS_LOG,
		BEAGLE_FLAG_VECTOR_NONE,
		BEAGLE_FLAG_VECTOR_SSE,
		BEAGLE_FLAG_THREADING_NONE,
		BEAGLE_FLAG_THREADING_OPENMP
	};
	int k = 0;
	for (int i = 0; i < 20; i++){
		if (inFlags & flags[i]){
			if (k%4 == 0 && k > 0)
				printf("\n                       ");
			printf(" %s", names[i]);
			k++;
		}
	}
}
#endif

double effectiveSampleSize(const std::vector<double>& logWeights)
{
    assert(!logWeights.empty());
    const double maxWeight = *std::max_element(logWeights.begin(), logWeights.end());
    std::vector<double> v(logWeights.size());

    auto normalize = [maxWeight](const double w) { return w - maxWeight; };
    std::transform(logWeights.begin(), logWeights.end(), v.begin(), normalize);

    auto dexp = [](const double d) { return std::exp(d); };
    std::transform(v.begin(), v.end(), v.begin(), dexp);

    const double s = std::accumulate(v.begin(), v.end(), 0.0);

    auto square = [maxWeight](const double w) { return std::pow(w, 2); };
    std::transform(v.begin(), v.end(), v.begin(), square);
    const double ssq = std::accumulate(v.begin(), v.end(), 0.0);

    return std::exp(-std::log(ssq) + 2 * std::log(s));
}
    
double logit(double x){
    return std::log(x/(1.0-x));
}
    
double logitinv(double x){
    return 1./(1.0+std::exp(-x));
}
} // namespace particle
} // namespace sts
