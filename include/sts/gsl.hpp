#ifndef STS_GSL_HPP
#define STS_GSL_HPP

#include <functional>
#include <cstdlib>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

namespace sts
{
namespace gsl
{
/// Convert a <c>std::function<double(double)></c> to a value usable with GSL.
/// \param d Input value
/// \param data a std::function<double(double)>
double std_func_to_gsl_function(double d, void *data)
{
    std::function<double(double)> *func = reinterpret_cast<std::function<double(double)>*>(data);
    return (*func)(d);
}

double minimize(const std::function<double(double)> fn,
                double m = 0.5,
                double a = 0,
                double b = 1,
                const int max_iter = 100,
                const gsl_min_fminimizer_type *min_type = gsl_min_fminimizer_brent)
{
    int iter = 0, status;
    gsl_min_fminimizer *s;
    gsl_function gsl_fn;

    gsl_fn.function = &std_func_to_gsl_function;
    gsl_fn.params = reinterpret_cast<void*>(&fn);

    s = gsl_min_fminimizer_alloc(min_type);
    gsl_min_fminimizer_set(s, &gsl_fn, m, a, b);

    do {
        iter++;
        status = gsl_min_fminimizer_iterate(s);

        m = gsl_min_fminimizer_x_minimum(s);
        a = gsl_min_fminimizer_x_lower(s);
        b = gsl_min_fminimizer_x_upper(s);

        status = gsl_min_test_interval(a, b, 0.001, 0.0);
    } while(status == GSL_CONTINUE && iter < max_iter);
    gsl_min_fminimizer_free(s);
    return m;
}

} // namespace gsl
} // namespace sts

#endif // STS_GSL_HPP
