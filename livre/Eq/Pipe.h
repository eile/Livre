
/* Copyright (c) 2006-2011, Stefan Eilemann <eile@equalizergraphics.com>
 *                    2011, Maxim Makhinya  <maxmah@gmail.com>
 *                    2012, David Steiner   <steiner@ifi.uzh.ch>
 *                    2013, Ahmet Bilgili   <ahmet.bilgili@epfl.ch>
 */

#ifndef _Pipe_h_
#define _Pipe_h_

#include <eq/pipe.h>
#include <livre/Eq/types.h>

namespace livre
{

namespace detail
{
class Pipe;
}

/**
 * The Pipe class is a standard EQ abstraction for OpenGL contexts.
 */
class Pipe : public eq::Pipe
{
public:

    /**
     * @param parent Parent node that owns the pipe.
     */
    Pipe( eq::Node* parent );
    ~Pipe();

    /**
     * @return The frame data.
     */
    ConstFrameDataPtr getFrameData( ) const;

private:

    void frameStart( const eq::uint128_t& frameId,
                     const uint32_t frameNumber ) final;
    bool configInit( const eq::uint128_t& initId ) final;
    bool configExit() final;

    detail::Pipe* _impl;
};

}

#endif //_Pipe_h_

