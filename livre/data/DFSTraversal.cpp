/* Copyright (c) 2011-2017, EPFL/Blue Brain Project
 *                     Ahmet Bilgili <ahmet.bilgili@epfl.ch>
 *
 * This file is part of Livre <https://github.com/BlueBrain/Livre>
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License version 3.0 as published
 * by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <livre/data/DFSTraversal.h>
#include <livre/data/DataSource.h>
#include <livre/data/NodeVisitor.h>

namespace livre
{
struct DFSTraversal::Impl
{
public:
    Impl(const DataSource& s)
        : source(s)
    {
    }

    void traverse(const NodeId& nodeId, const uint32_t depth,
                  livre::NodeVisitor& visitor)
    {
        assert(depth > 0);

        if (!visitor.visit(source.getNode(nodeId)))
            return;

        if (depth == 1)
            return;

        const NodeIds& nodeIds = nodeId.getChildren();
        for (const NodeId& childNodeId : nodeIds)
            traverse(childNodeId, depth - 1, visitor);
    }

    const DataSource& source;
};

DFSTraversal::DFSTraversal(const DataSource& source)
    : _impl(new DFSTraversal::Impl(source))
{
}

DFSTraversal::~DFSTraversal() = default;

void DFSTraversal::traverse(NodeVisitor& visitor, const uint32_t timeStep)
{
    const auto& rootNode = _impl->source.getVolumeInfo().rootNode;
    const auto depth = rootNode.getDepth();
    const Vector3ui& blockSize = rootNode.getBlockSize();

    visitor.visitPre();
    for (uint32_t x = 0; x < blockSize.x(); ++x)
        for (uint32_t y = 0; y < blockSize.y(); ++y)
            for (uint32_t z = 0; z < blockSize.z(); ++z)
            {
                _impl->traverse(NodeId(0, Vector3ui(x, y, z), timeStep), depth,
                                visitor);
            }
    visitor.visitPost();
}
}
