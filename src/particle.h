#ifndef STS_PARTICLE_PARTICLE_H
#define STS_PARTICLE_PARTICLE_H

#include <memory>

namespace sts
{
namespace particle
{
class State;
typedef std::shared_ptr<State> Particle;

} // namespace particle
} // namespace sts

#endif // STS_PARTICLE_PARTICLE_H

