/*
 * cube.h
 *
 *  Created on: Dec 4, 2011
 *      Author: Jurgen
 */

#ifndef MESH_H
#define FLAG_H


#include "err.h"
#include "entity.h"
#include "brush.h"
#include "texture.h"
#include "entitypool.h"
#include "allocator.h"

#include <vector>

class Flag : public Entity
{
public:
    enum {
        BLEND_COLOR_F = 1<<RenderState::USER_B
    };
    enum {
        TEXTURE,
        NORMAL,
        SPECULAR
    };

private:
    GLuint      m_VboID;
protected:
    std::vector< BrushPtr > m_Assets;    // all mmaped data. Don't worry about pixel storage memory here!
    TexturePtr  m_Texture;
    TexturePtr  m_NormalMap;
    TexturePtr  m_SpecularNap;

    typedef std::vector<Vector, Allocator<Vector>> VertexVector;

    VertexVector m_VertexBuffer;      // linear vertex buffer - custom allocator
    VertexVector m_TexCoordBuffer;    // linear texture coord buffer - custom allocator (4 component...for now)
    std::vector<int>  m_IndexArray;   // standard array to map vertices to tris

    boost::shared_ptr<MemoryPool> m_MemoryPool;

public:
    Flag( const std::vector< BrushPtr >& assets );

    virtual ~Flag();

protected:
    virtual bool DoInitialize( Renderer* renderer ) throw(std::exception);

    virtual void DoRender() throw(std::exception);

    virtual void DoUpdate( float ticks ) throw(std::exception);

};

#if 0
#include <render/drawable.h>

#include <vector>

namespace wave
{
	class WaveMesh : public pei::Drawable
	{
	    std::vector<GLfloat> texCoords;
	    double 	 amp;
	    float 	 xrot, yrot, zrot;
	    double 	 m_TimeEllapsed;
	    double	 m_Speed;
	public:
		WaveMesh();

		virtual ~WaveMesh();

		virtual pei::SurfacePtr SetSurface( pei::SurfacePtr s, int n = 0 );

	protected:
		virtual void OnUpdateAnimation( const pei::RenderProfilePtr& profile, double time );

        virtual void OnUpdateState( const pei::RenderProfilePtr& profile );

		virtual void OnDraw( const pei::RenderProfilePtr& profile, const pei::SurfacePtr& buffer, const pei::RenderParam& param  );
	};

}

#endif

#endif /* MESH_H */
