#ifndef OCTOMAP_LABEL_COUNTING_OCTREE_HH
#define OCTOMAP_LABEL_COUNTING_OCTREE_HH

#include <stdio.h>
#include <valarray>
#include "OcTreeBase.h"
#include "OcTreeDataNode.h"

namespace octomap
{
  /**
   * An Octree-node which stores an internal counter per node / volume.
   *
   * Count is recursive, parent nodes have the summed count of their
   * children.
   *
   * \note In our mapping system this data structure is used in
   *       LabelCountingOcTree in the sensor model only
   */
  class LabelCountingOcTreeNode : public OcTreeDataNode<std::valarray<unsigned int> >
  {
  public:
    LabelCountingOcTreeNode();
    ~LabelCountingOcTreeNode();

    inline std::valarray<unsigned int> getCount() const { return getValue(); }
    inline void decreaseCount(unsigned int label, unsigned int n_label)
    {
      if (!is_value_initialized_)
      {
        value = std::valarray<unsigned int>(static_cast<unsigned int>(0), n_label);
        is_value_initialized_ = true;
      }
      if (0 <= label && label < value.size())
      {
        if (value[label] > 0)
        {
          value[label]--;
        }
      }
      else
      {
        printf("label value must satisfy: 0 <= label < %zu\n", value.size());
      }
    }
    inline void increaseCount(unsigned int label, unsigned int n_label)
    {
      if (!is_value_initialized_)
      {
        value = std::valarray<unsigned int>(static_cast<unsigned int>(0), n_label);
        is_value_initialized_ = true;
      }
      if (0 <= label && label < value.size())
      {
        value[label]++;
      }
      else
      {
        printf("label value must satisfy: 0 <= label < %zu\n", value.size());
      }
    }
    inline void setCount(std::valarray<unsigned int> c) { this->setValue(c); }

  protected:
    bool is_value_initialized_;
    unsigned int n_label_;
  };

  /**
   * An AbstractOcTree which stores an internal counter per node / volume.
   *
   * Count is recursive, parent nodes have the summed count of their
   * children.
   *
   * \note Was only used internally, not used anymore
   */
  class LabelCountingOcTree : public OcTreeBase <LabelCountingOcTreeNode> {

  public:
    /// Default constructor, sets resolution of leafs
    LabelCountingOcTree(double resolution, unsigned int n_label);
    virtual LabelCountingOcTreeNode* updateNode(
      const point3d& value, int label, const bool hit=true, const bool reset=false);
    LabelCountingOcTreeNode* updateNode(
      const OcTreeKey& k, int label, const bool hit=true, const bool reset=false);
    void getCentersMinHits(point3d_list& node_centers, std::vector<unsigned int>& labels, unsigned int min_hits) const;

  protected:

    // number of labels and size of the valarray (value)
    unsigned int n_label_;

    void getCentersMinHitsRecurs(point3d_list& node_centers,
                                 std::vector<unsigned int>& labels,
                                 unsigned int min_hits,
                                 unsigned int max_depth,
                                 LabelCountingOcTreeNode* node, unsigned int depth,
                                 const OcTreeKey& parent_key) const;

    /**
     * Static member object which ensures that this OcTree's prototype
     * ends up in the classIDMapping only once. You need this as a
     * static member in any derived octree class in order to read .ot
     * files through the AbstractOcTree factory. You should also call
     * ensureLinking() once from the constructor.
     */
    class StaticMemberInitializer{
       public:
         StaticMemberInitializer() {
           LabelCountingOcTree* tree = new LabelCountingOcTree(0.1, 1);
           tree->clearKeyRays();
           AbstractOcTree::registerTreeType(tree);
         }

         /**
         * Dummy function to ensure that MSVC does not drop the
         * StaticMemberInitializer, causing this tree failing to register.
         * Needs to be called from the constructor of this octree.
         */
         void ensureLinking() {};
    };
    /// static member to ensure static initialization (only once)
    static StaticMemberInitializer labelCountingOcTreeMemberInit;
  };

}  // namespace octomap


#endif  // OCTOMAP_LABEL_COUNTING_OCTREE_HH
