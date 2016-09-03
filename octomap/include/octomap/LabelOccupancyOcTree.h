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

#ifndef OCTOMAP_LABELOCCUPANCYOCTREE_H
#define OCTOMAP_LABELOCCUPANCYOCTREE_H

#include "OccupancyOcTreeBase.h"
#include "octomap_types.h"
#include "octomap_utils.h"
#include "OcTreeDataNode.h"
#include <limits>
#include <valarray>

namespace octomap
{

  class LabelOccupancyOcTree : public OccupancyOcTreeBase<LabelOccupancyOcTreeNode>
  {
  public:
    /// Default constructor, sets resolution of leafs and number of labels
    LabelOccupancyOcTree(double resolution, int n_label);
    virtual ~LabelOccupancyOcTree() {};

    /// virtual constructor: creates a new object of same type
    /// (Covariant return type requires an up-to-date compiler)
    LabelOccupancyOcTree* create() const { return new LabelOccupancyOcTree(resolution, n_label); }

    std::string getTreeType() const { return "LabelOccupancyOcTree"; }

  protected:
    // Number of labels
    int n_label;

    /**
     * Static member object which ensures that this OcTree's prototype
     * ends up in the classIDMapping only once. You need this as a
     * static member in any derived octree class in order to read .ot
     * files through the AbstractOcTree factory. You should also call
     * ensureLinking() once from the constructor.
     */
    class StaticMemberInitializer
    {
    public:
      StaticMemberInitializer()
      {
        LabelOccupancyOcTree* tree = new LabelOccupancyOcTree(0.1, 1);
        tree->clearKeyRays();
        AbstractOcTree::registerTreeType(tree);
      }

      /**
       * Dummy function to ensure that MSVC does not drop the
       * StaticMemberInitializer, causing this tree failing to register.
       * Needs to be called from the constructor of this octree.
       */
      void ensureLinking() {}
    };

    /// to ensure static initialization (only once)
    static StaticMemberInitializer ocTreeMemberInit;
  };

} // namespace octomap

#endif  // OCTOMAP_LABELOCCUPANCYOCTREE_H
