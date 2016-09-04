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

#ifndef OCTOMAP_LABELOCCUPANCYOCTREENODE_H
#define OCTOMAP_LABELOCCUPANCYOCTREENODE_H

#include "octomap_types.h"
#include "octomap_utils.h"
#include "OcTreeDataNode.h"
#include <limits>
#include <valarray>

namespace octomap
{

  /**
   * Nodes to be used in OcTree. They represent 3d occupancy grid cells.
   * "value" stores their log-odds occupancy.
   *
   * Note: If you derive a class (directly or indirectly) from OcTreeNode or
   * OcTreeDataNode, you have to implement (at least) the following functions:
   * createChild(), getChild(), getChild() const, expandNode() to avoid slicing
   * errors and memory-related bugs.
   * See ColorOcTreeNode in ColorOcTree.h for an example.
   *
   */
  class LabelOccupancyOcTreeNode : public OcTreeDataNode<std::valarray<double> >
  {
  public:
    LabelOccupancyOcTreeNode();
    ~LabelOccupancyOcTreeNode();

    // -- node occupancy  ----------------------------

    /// \return occupancy probability of node
    inline std::valarray<double> getOccupancy() const
    {
      std::valarray<double> occupancy(value.size());
      for (size_t i=0; i < value.size(); i++)
      {
        occupancy[i] = probability(value[i]);
      }
      return occupancy;
    }

    /// \return log odds representation of occupancy probability of node
    inline std::valarray<double> getLogOdds() const { return value; }
    /// sets log odds occupancy of node
    inline void setLogOdds(std::valarray<double> l) { value = l; }

    /**
     * @return mean of all children's occupancy probabilities, in log odds
     */
    std::valarray<double> getMeanChildLogOdds() const;

    /**
     * @return max of all children's occupancy probabilities, in log odds
     */
    std::valarray<double> getMaxChildLogOdds() const;

    /// update this node's occupancy according to its children's maximum occupancy
    inline void updateOccupancyChildren()
    {
      this->setLogOdds(this->getMeanChildLogOdds());  // conservative
    }

    /// adds p to the node's logOdds value (with no boundary / threshold checking!)
    void addValue(const std::valarray<double> p);


  protected:
    // "value" stores log odds occupancy probability
  };

} // namespace octomap

#endif  // OCTOMAP_LABELOCCUPANCYOCTREENODE_H
