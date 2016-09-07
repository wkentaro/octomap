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

#include <octomap/LabelOccupancyOcTreeNode.h>

namespace octomap
{

  LabelOccupancyOcTreeNode::LabelOccupancyOcTreeNode()
    : OcTreeDataNode<std::valarray<float> >(),
      is_value_initialized(false)
  {
  }

  LabelOccupancyOcTreeNode::~LabelOccupancyOcTreeNode()
  {
  }

  std::valarray<float> LabelOccupancyOcTreeNode::getMeanChildLogOdds() const
  {
    std::valarray<float> mean;

    uint8_t c = 0;
    if (children != NULL)
    {
      for (unsigned int i=0; i < 8; i++)
      {
        if (children[i] != NULL)
        {
          std::valarray<float> occupancy = static_cast<LabelOccupancyOcTreeNode*>(children[i])->getOccupancy();
          if (mean.size() == 0)
          {
            mean = std::valarray<float>(occupancy.size());
          }
          else
          {
            mean += occupancy;
          }
        }
      }
    }

    if (c > 0)
    {
      mean /= static_cast<float>(c);
    }

    return log(mean / (std::valarray<float>(1, mean.size()) - mean));
  }

  std::valarray<float> LabelOccupancyOcTreeNode::getMaxChildLogOdds() const
  {
    std::valarray<float> maxLogOdds;
    float maxOccupancy = -std::numeric_limits<float>::max();

    if (children !=NULL){
      for (unsigned int i=0; i<8; i++) {
        if (children[i] != NULL) {
          std::valarray<float> l = static_cast<LabelOccupancyOcTreeNode*>(children[i])->getLogOdds(); // TODO check if works generally
          if (l.max() > maxOccupancy)
          {
            maxOccupancy = l.max();
            maxLogOdds = l;
          }
        }
      }
    }
    return maxLogOdds;
  }

  void LabelOccupancyOcTreeNode::addValue(const std::valarray<float>& logOdds)
  {
    if (!is_value_initialized)
    {
      value = std::valarray<float>(0.0, logOdds.size());
      is_value_initialized = true;
    }
    assert(value.size() == logOdds.size());
    value += logOdds;
  }

}  // namespace octomap
