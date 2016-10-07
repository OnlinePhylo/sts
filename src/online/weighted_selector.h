#ifndef STS_ONLINE_WEIGHTED_SELECTOR_H
#define STS_ONLINE_WEIGHTED_SELECTOR_H

#include <smctc.hh>
#include <stdexcept>
#include <vector>
#include <numeric>

namespace sts { namespace online {

template<typename T>
class WeightedSelector
{
public:
    WeightedSelector(smc::rng& rng): _rng(rng){};
    WeightedSelector(smc::rng& rng, const std::vector<T>& v, const std::vector<double>& weights);

    void push_back(T t, const double weight);

    const T& choice() const;
    T& choice();

    inline size_t size() const { return values.size(); };
private:
    std::vector<T> values;
    std::vector<double> weights;
    smc::rng& _rng;

    size_t choose_index() const;
};

template<typename T>
WeightedSelector<T>::WeightedSelector(smc::rng& rng, const std::vector<T>& v,
                                  const std::vector<double>& weights) :
    values(v),
    weights(weights),
    _rng(rng)
{
    if(values.size() != weights.size())
        throw std::runtime_error("Number of values does not match number of weights.");
}


template<typename T>
void WeightedSelector<T>::push_back(T t, const double weight)
{
    values.push_back(t);
    weights.push_back(weight);
}

template<typename T>
size_t WeightedSelector<T>::choose_index() const
{
    std::vector<unsigned int> indices(size());

    _rng.Multinomial(1, size(), weights.data(), indices.data());
    size_t i = 0;
    for(; i < size(); i++) {
        if(indices[i] > 0) {
            assert(indices[i] == 1);
            break;
        }
    }
    return i;
}


template<typename T>
const T& WeightedSelector<T>::choice() const
{
    return values[choose_index()];
}

template<typename T>
T& WeightedSelector<T>::choice()
{
    return values[choose_index()];
}


}} // Namespaces

#endif
