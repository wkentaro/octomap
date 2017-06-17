#include <cassert>
#include <octomap/LabelCountingOcTree.h>

namespace octomap
{
  /// implementation of LabelCountingOcTreeNode  ----------------------------------

  LabelCountingOcTreeNode::LabelCountingOcTreeNode()
    : OcTreeDataNode<std::valarray<unsigned int> >(), is_value_initialized_(false)
  {
  }

  LabelCountingOcTreeNode::~LabelCountingOcTreeNode()
  {
  }

  /// implementation of LabelCountingOcTree  --------------------------------------

  LabelCountingOcTree::LabelCountingOcTree(double resolution, unsigned int n_label)
   : OcTreeBase<LabelCountingOcTreeNode>(resolution), n_label_(n_label)
  {
    labelCountingOcTreeMemberInit.ensureLinking();
  }

  LabelCountingOcTreeNode* LabelCountingOcTree::updateNode(const point3d& value, unsigned int label)
  {
    OcTreeKey key;
    if (!coordToKeyChecked(value, key)) return NULL;
    return updateNode(key, label);
  }

  // Note: do not inline this method, will decrease speed (KMW)
  LabelCountingOcTreeNode* LabelCountingOcTree::updateNode(const OcTreeKey& k, unsigned int label)
  {
    if (root == NULL)
    {
      root = new LabelCountingOcTreeNode();
      tree_size++;
    }
    LabelCountingOcTreeNode* curNode(root);
    curNode->increaseCount(label, n_label_);

    // follow or construct nodes down to last level...
    for (int i=(tree_depth-1); i>=0; i--)
    {
      unsigned int pos = computeChildIdx(k, i);

      // requested node does not exist
      if (!nodeChildExists(curNode, pos))
      {
        createNodeChild(curNode, pos);
      }
      // descent tree
      curNode = getNodeChild(curNode, pos);
      curNode->increaseCount(label, n_label_); // modify traversed nodes
    }
    return curNode;
  }

  void LabelCountingOcTree::getCentersMinHits(point3d_list& node_centers,
                                              std::vector<unsigned int>& labels,
                                              unsigned int min_hits) const
  {
    OcTreeKey root_key;
    root_key[0] = root_key[1] = root_key[2] = this->tree_max_val;
    getCentersMinHitsRecurs(node_centers, labels, min_hits, this->tree_depth, this->root, 0, root_key);
  }

  void LabelCountingOcTree::getCentersMinHitsRecurs(point3d_list& node_centers,
                                                    std::vector<unsigned int>& labels,
                                                    unsigned int min_hits,
                                                    unsigned int max_depth,
                                                    LabelCountingOcTreeNode* node, unsigned int depth,
                                                    const OcTreeKey& parent_key) const
  {
    if (depth < max_depth && nodeHasChildren(node))
    {
      key_type center_offset_key = this->tree_max_val >> (depth + 1);
      OcTreeKey search_key;

      for (unsigned int i = 0; i < 8; ++i)
      {
        if (nodeChildExists(node, i))
        {
          computeChildKey(i, center_offset_key, parent_key, search_key);
          getCentersMinHitsRecurs(node_centers, labels, min_hits, max_depth,
                                  getNodeChild(node, i), depth + 1, search_key);
        }
      }
    }
    else
    {
      std::valarray<unsigned int> count = node->getCount();
      unsigned int max_hits = 0;
      unsigned int max_hits_label = 0;
      for (size_t label = 0; label < count.size(); label++)
      {
        if (count[label] > max_hits) {
          max_hits = count[label];
          max_hits_label = label;
        }
      }
      // max level reached
      if (max_hits >= min_hits)
      {
        node_centers.push_back(this->keyToCoord(parent_key, depth));
        labels.push_back(max_hits_label);
      }
    }
  }

  LabelCountingOcTree::StaticMemberInitializer LabelCountingOcTree::labelCountingOcTreeMemberInit;

} // namespace octomap
