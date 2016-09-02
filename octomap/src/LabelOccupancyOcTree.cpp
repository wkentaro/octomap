/*
 * OctoMap - An Efficient Probabilistic 3D Mapping Framework Based on Octrees
 * http://octomap.github.com/
 *
 * Copyright (c) 2016, Kentaro Wada.
 * All rights reserved.
 * License: New BSD
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Freiburg nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <algorithm>
#include <valarray>
#include <bitset>
#include <cassert>
#include <math.h>
#include <fstream>
#include <utility>
#include <stdlib.h>
#include <inttypes.h>

#include "octomap/LabelOccupancyOcTreeNode.h"
#include "octomap/LabelOccupancyOcTree.h"

namespace octomap
{

  LabelOccupancyOcTree::LabelOccupancyOcTree(double resolution, int n_label)
    : public OccupancyOcTreeBase<LabelOccupancyOcTreeNode>(resolution)
  {
    this->n_label = n_label;
    ocTreeMemberInit.ensureLinking();
  };

  OcTree::StaticMemberInitializer OcTree::ocTreeMemberInit;

  // Overwrites methods of OccupancyOcTreeBase for using multiple label occupancies

  LabelOccupancyOcTreeNode* LabelOccupancyOcTree::setNodeValue(
      const OcTreeKey& key, std::valarray<double> log_odds_value, bool lazy_eval)
  {
    assert(this->n_label == log_odds_value.size());

    // clamp log odds within range:
    for (size_t i=0; i < log_odds_value.size(); i++)
    {
      log_odds_value[i] = std::min(std::max(log_odds_value[i], this->clamping_thres_min), this->clamping_thres_max);
    }

    bool createdRoot = false;
    if (this->root == NULL)
    {
      this->root = new LabelOccupancyOcTreeNode();
      this->tree_size++;
      createdRoot = true;
    }

    return setNodeValueRecurs(this->root, createdRoot, key, 0, log_odds_value, lazy_eval);
  }

  LabelOccupancyOcTreeNode* LabelOccupancyOcTree::setNodeValue(
      const point3d& value, std::valarray<double> log_odds_value, bool lazy_eval)
  {
    assert(this->n_label == log_odds_value.size());

    OcTreeKey key;
    if (!this->coordToKeyChecked(value, key))
    {
      return NULL;
    }

    return setNodeValue(key, log_odds_value, lazy_eval);
  }

  LabelOccupancyOcTreeNode* LabelOccupancyOcTree::setNodeValue(
      double x, double y, double z, std::valarray<double> log_odds_value, bool lazy_eval)
  {
    assert(this->n_label == log_odds_value.size());

    OcTreeKey key;
    if (!this->coordToKeyChecked(x, y, z, key))
    {
      return NULL;
    }

    return setNodeValue(key, log_odds_value, lazy_eval);
  }

  LabelOccupancyOcTreeNode* LabelOccupancyOcTree::updateNode(
      const OcTreeKey& key, std::valarray<double> log_odds_update, bool lazy_eval)
  {
    assert(this->n_label == log_odds_update.size());

    bool createdRoot = false;
    if (this->root == NULL)
    {
      this->root = new LabelOccupancyOcTree();
      this->tree_size++;
      createdRoot = true;
    }

    return updateNodeRecurs(this->root, createdRoot, key, 0, log_odds_update, lazy_eval);
  }

  LabelOccupancyOcTreeNode* LabelOccupancyOcTree::updateNode(
      const point3d& value, std::valarray<double> log_odds_update, bool lazy_eval)
  {
    assert(this->n_label == log_odds_update.size());

    OcTreeKey key;
    if (!this->coordToKeyChecked(value, key))
    {
      return NULL;
    }

    return updateNode(key, log_odds_update, lazy_eval);
  }

  LabelOccupancyOcTreeNode* LabelOccupancyOcTree::updateNode(
      double x, double y, double z, std::valarray<double> log_odds_update, bool lazy_eval)
  {
    assert(this->n_label == log_odds_update.size());

    OcTreeKey key;
    if (!this->coordToKeyChecked(x, y, z, key))
    {
      return NULL;
    }

    return updateNode(key, log_odds_update, lazy_eval);
  }

  LabelOccupancyOcTreeNode* LabelOccupancyOcTree::updateNode(const OcTreeKey& key, bool occupied, bool lazy_eval)
  {
    std::valarray<double> logOdds;
    if (occupied)
    {
      logOdds = std::valarray<double>(this->prob_hit_log, this->n_label);
    }
    else
    {
      logodds = std::valarray<double>(this->prob_miss_log, this->n_label);
    }

    return updateNode(key, logOdds, lazy_eval);
  }

  LabelOccupancyOcTreeNode* LabelOccupancyOcTree::updateNode(const point3d& value, bool occupied, bool lazy_eval)
  {
    OcTreeKey key;
    if (!this->coordToKeyChecked(value, key))
    {
      return NULL;
    }

    return updateNode(key, occupied, lazy_eval);
  }

  LabelOccupancyOcTreeNode* LabelOccupancyOcTree::updateNode(
      double x, double y, double z, bool occupied, bool lazy_eval)
  {
    OcTreeKey key;
    if (!this->coordToKeyChecked(x, y, z, key))
    {
      return NULL;
    }

    return updateNode(key, occupied, lazy_eval);
  }

  inline bool LabelOccupancyOcTree::isNodeOccupied(const LabelOccupancyOcTreeNode* occupancyNode) const
  {
    return (occupancyNode->getLogOdds().max() >= this->occ_prob_thres_log);
  }

  inline bool LabelOccupancyOcTree::isNodeOccupied(const LabelOccupancyOcTreeNode& occupancyNode) const
  {
    return (occupancyNode.getLogOdds().max() >= this->occ_prob_thres_log);
  }

  inline bool LabelOccupancyOcTree::isNodeAtThreshold(const LabelOccupancyOcTreeNode* occupancyNode) const
  {
    throw std::runtime_error("LabelOccupancyOcTree::isNodeAtThreshold is not supported.");
  }

  inline bool LabelOccupancyOcTree::isNodeAtThreshold(const OcTreeNode& occupancyNode) const
  {
    throw std::runtime_error("LabelOccupancyOcTree::isNodeAtThreshold is not supported.");
  }

  LabelOccupancyOcTreeNode* LabelOccupancyOcTree::updateNodeRecurs(
      LabelOccupancyOcTreeNode* node, bool node_just_created, const OcTreeKey& key,
      unsigned int depth, const std::valarray<double> log_odds_update, bool lazy_eval)
  {
    bool created_node = false;

    assert(node);

    // follow down to last level
    if (depth < this->tree_depth)
    {
      unsigned int pos = computeChildIdx(key, this->tree_depth - 1 - depth);
      if (!this->nodeChildExists(node, pos))
      {
        // child does not exist, but maybe it's a pruned node?
        if (!this->nodeHasChildren(node) && !node_just_created)
        {
          // current node does not have children AND it is not a new node
          // -> expand pruned node
          this->expandNode(node);
        }
        else
        {
          // not a pruned node, create requested child
          this->createNodeChild(node, pos);
          created_node = true;
        }
      }

      if (lazy_eval)
      {
        return updateNodeRecurs(this->getNodeChild(node, pos), created_node, key, depth+1, log_odds_update, lazy_eval);
      }
      else
      {
        LabelOccupancyOcTreeNode* retval = updateNodeRecurs(
            this->getNodeChild(node, pos), created_node, key, depth+1, log_odds_update, lazy_eval);
        // prune node if possible, otherwise set own probability
        // note: combining both did not lead to a speedup!
        if (this->pruneNode(node))
        {
          // return pointer to current parent (pruned), the just updated node no longer exists
          retval = node;
        }
        else
        {
          node->updateOccupancyChildren();
        }

        return retval;
      }
    }
    // at last level, update node, end of recursion
    else
    {
      if (use_change_detection)
      {
        bool occBefore = this->isNodeOccupied(node);
        updateNodeLogOdds(node, log_odds_update);

        if (node_just_created)  // new node
        {
          changed_keys.insert(std::pair<OcTreeKey, bool>(key, true));
        }
        else if (occBefore != this->isNodeOccupied(node))  // occupancy changed, track it
        {
          KeyBoolMap::iterator it = changed_keys.find(key);
          if (it == changed_keys.end())
          {
            changed_keys.insert(std::pair<OcTreeKey, bool>(key, false));
          }
          else if (it->second == false)
          {
            changed_keys.erase(it);
          }
        }
      }
      else
      {
        updateNodeLogOdds(node, log_odds_update);
      }
      return node;
    }
  }

  LabelOccupancyOcTreeNode* LabelOccupancyOcTree::setNodeValueRecurs(
      LabelOccupancyOcTreeNode* node, bool node_just_created, const OcTreeKey& key,
      unsigned int depth, const std::valarray<double> log_odds_value, bool lazy_eval)
  {
    bool created_node = false;

    assert(node);

    // follow down to last level
    if (depth < this->tree_depth)
    {
      unsigned int pos = computeChildIdx(key, this->tree_depth -1 - depth);
      if (!this->nodeChildExists(node, pos))
      {
        // child does not exist, but maybe it's a pruned node?
        if (!this->nodeHasChildren(node) && !node_just_created)
        {
          // current node does not have children AND it is not a new node
          // -> expand pruned node
          this->expandNode(node);
        }
        else
        {
          // not a pruned node, create requested child
          this->createNodeChild(node, pos);
          created_node = true;
        }
      }

      if (lazy_eval)
      {
        return setNodeValueRecurs(this->getNodeChild(node, pos), created_node, key,
                                  depth+1, log_odds_value, lazy_eval);
      }
      else
      {
        LabelOccupancyOcTreeNode* retval = setNodeValueRecurs(
            this->getNodeChild(node, pos), created_node, key, depth+1, log_odds_value, lazy_eval);
        // prune node if possible, otherwise set own probability
        // note: combining both did not lead to a speedup!
        if (this->pruneNode(node))
        {
          // return pointer to current parent (pruned), the just updated node no longer exists
          retval = node;
        }
        else
        {
          node->updateOccupancyChildren();
        }

        return retval;
      }
    }
    // at last level, update node, end of recursion
    else
    {
      if (use_change_detection)
      {
        bool occBefore = this->isNodeOccupied(node);
        node->setLogOdds(log_odds_value);

        if (node_just_created)  // new node
        {
          changed_keys.insert(std::pair<OcTreeKey, bool>(key, true));
        }
        else if (occBefore != this->isNodeOccupied(node))  // occupancy changed, track it
        {
          KeyBoolMap::iterator it = changed_keys.find(key);
          if (it == changed_keys.end())
          {
            changed_keys.insert(std::pair<OcTreeKey, bool>(key, false));
          }
          else if (it->second == false)
          {
            changed_keys.erase(it);
          }
        }
      }
      else
      {
        node->setLogOdds(log_odds_value);
      }
      return node;
    }
  }

  void LabelOccupancyOcTree::updateNodeLogOdds(
      LabelOccupancyOcTreeNode* occupancyNode, const std::valarray<double> update) const
  {
    occupancyNode->addValue(update);
    std::valarray<double> logOdds = occupancyNode->getLogOdds();
    for (size_t i=0; i < logOdds.size(); i++)
    {
      if (logOdds[i] < this->clamping_thres_min)
      {
        logOdds[i] = this->clamping_thres_min;
      }
      if (logOdds[i] > this->clamping_thres_max)
      {
        logOdds[i] = this->clamping_thres_max;
      }
    }
  }

  void LabelOccupancyOcTree::integrateHit(LabelOccupancyOcTreeNode* occupancyNode) const
  {
    std::valarray<double> logOdds = std::valarray(this->prob_hit_log, this->n_label);
    updateNodeLogOdds(occupancyNode, logOdds);
  }

  void LabelOccupancyOcTree::integrateMiss(LabelOccupancyOcTreeNode* occupancyNode) const
  {
    std::valarray<double> logOdds = std::valarray(this->prob_miss_log, this->n_label);
    updateNodeLogOdds(occupancyNode, logOdds);
  }

  void LabelOccupancyOcTree::nodeToMaxLikelihood(LabelOccupancyOcTreeNode* occupancyNode) const
  {
    std::valarray<double> logOdds;
    if (this->isNodeOccupied(occupancyNode))
    {
      logOdds = std::valarray(this->prob_hit_log, this->n_label);
    }
    else
    {
      logOdds = std::valarray(this->prob_miss_log, this->n_label);
    }
    occupancyNode->setLogOdds(logOdds);
  }

  void LabelOccupancyOcTree::nodeToMaxLikelihood(LabelOccupancyOcTreeNode& occupancyNode) const
  {
    std::valarray<double> logOdds;
    if (this->isNodeOccupied(occupancyNode))
    {
      logOdds = std::valarray(this->prob_hit_log, this->n_label);
    }
    else
    {
      logOdds = std::valarray(this->prob_miss_log, this->n_label);
    }
    occupancyNode.setLogOdds(logOdds);
  }

}  // namespace octomap
