#include "online_util.h"

#include <gsl/gsl_statistics_double.h>
namespace sts { namespace online {
SummaryStatistics summarize(std::vector<double> values)
{
    std::sort(values.begin(), values.end()); // For median
    const double* arr = values.data();
    const size_t n = values.size();
    SummaryStatistics result;
    result.mean = gsl_stats_mean(arr, 1, n);
    result.median = gsl_stats_median_from_sorted_data(arr, 1, n);
    result.lower95 = gsl_stats_quantile_from_sorted_data(arr, 1, n, 0.025);
    result.upper95 = gsl_stats_quantile_from_sorted_data(arr, 1, n, 0.975);
    result.sd = gsl_stats_sd(arr, 1, n);
    result.min = values.front();
    result.max = values.back();
    return result;
}
}}

namespace std {
std::ostream& operator<<(std::ostream& os, const sts::online::SummaryStatistics& value)
{
        return os << "min:    " << value.min << '\n'
                  << "2.5%:   " << value.lower95 << '\n'
                  << "median: " << value.median << '\n'
                  << "mean:   " << value.mean << '\n'
                  << "97.5%:  " << value.upper95 << '\n'
                  << "max:    " << value.max << '\n';
}
}


