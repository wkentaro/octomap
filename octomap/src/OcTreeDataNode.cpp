#include <valarray>

#include <octomap/OcTreeDataNode.h>

namespace octomap
{
  template <>
  bool OcTreeDataNode<std::valarray<unsigned int> >::operator==(
    const OcTreeDataNode<std::valarray<unsigned int> >& rhs) const
  {
    for (size_t i = 0; i < rhs.value.size(); i++)
    {
      if (rhs.value[i] != value[i])
      {
        return false;
      }
    }
    return true;
  }

} // namespace octomap
