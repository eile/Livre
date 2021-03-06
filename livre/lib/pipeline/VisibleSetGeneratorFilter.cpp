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

#include <livre/lib/configuration/VolumeRendererParameters.h>
#include <livre/lib/pipeline/VisibleSetGeneratorFilter.h>

#include <livre/core/cache/Cache.h>
#include <livre/core/pipeline/InputPort.h>
#include <livre/core/pipeline/PortData.h>
#include <livre/core/pipeline/Workers.h>
#include <livre/data/DFSTraversal.h>
#include <livre/data/DataSource.h>
#include <livre/data/SelectVisibles.h>

namespace livre
{
struct VisibleSetGeneratorFilter::Impl
{
    explicit Impl(const DataSource& dataSource)
        : _dataSource(dataSource)
    {
    }

    void execute(const FutureMap& input, PromiseMap& output) const
    {
        const UniqueFutureMap uniqueInputs(input.getFutures());

        const auto& frustum = uniqueInputs.get<Frustum>("Frustum");
        const auto& frame = uniqueInputs.get<uint32_t>("Frame");
        const auto& range = uniqueInputs.get<Range>("DataRange");
        const auto& params =
            uniqueInputs.get<VolumeRendererParameters>("Params");
        const auto& vp = uniqueInputs.get<PixelViewport>("Viewport");
        const auto& clipPlanes = uniqueInputs.get<ClipPlanes>("ClipPlanes");

        const uint32_t windowHeight = vp[3];
        const float sse = params.getScreenSpaceError();
        const uint32_t minLOD = params.getMinLod();
        const uint32_t maxLOD =
            std::min(params.getMaxLod(),
                     _dataSource.getVolumeInfo().rootNode.getDepth() - 1);

        SelectVisibles visitor(frustum, windowHeight, sse, minLOD, maxLOD,
                               range, clipPlanes);

        DFSTraversal traverser(_dataSource);
        traverser.traverse(visitor, frame);

        output.set("VisibleNodes", visitor.getVisibles());
        output.set("Params", params);
    }

    DataInfos getInputDataInfos() const
    {
        return {{"Frustum", getType<Frustum>()},
                {"Frame", getType<uint32_t>()},
                {"DataRange", getType<Range>()},
                {"Params", getType<VolumeRendererParameters>()},
                {"Viewport", getType<PixelViewport>()},
                {"ClipPlanes", getType<ClipPlanes>()}};
    }

    DataInfos getOutputDataInfos() const
    {
        return {{"VisibleNodes", getType<NodeIds>()},
                {"Params", getType<VolumeRendererParameters>()}};
    }

    const DataSource& _dataSource;
};

VisibleSetGeneratorFilter::VisibleSetGeneratorFilter(
    const DataSource& dataSource)
    : _impl(new VisibleSetGeneratorFilter::Impl(dataSource))
{
}

VisibleSetGeneratorFilter::~VisibleSetGeneratorFilter()
{
}

void VisibleSetGeneratorFilter::execute(const FutureMap& input,
                                        PromiseMap& output) const
{
    _impl->execute(input, output);
}

DataInfos VisibleSetGeneratorFilter::getInputDataInfos() const
{
    return _impl->getInputDataInfos();
}

DataInfos VisibleSetGeneratorFilter::getOutputDataInfos() const
{
    return _impl->getOutputDataInfos();
}
}
