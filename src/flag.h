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

class Surface : public Entity
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

    std::vector< BrushPtr > m_Assets;    // all mmaped data. Don't worry about pixel storage memory here!
    TexturePtr  m_Texture;
    TexturePtr  m_NormalMap;
    TexturePtr  m_SpecularNap;

    typedef std::vector<Vector, Allocator<Vector>> VertexVector;

    boost::shared_ptr<MemoryPool> m_MemoryPool;

    VertexVector m_VertexBuffer;      // linear vertex buffer - custom allocator
    VertexVector m_TexCoordBuffer;    // linear texture coord buffer - custom allocator (4 component...for now)
    std::vector<int>  m_IndexArray;   // standard array to map vertices to tris

    float       m_TimeEllapsed;
    float       m_Speed;

    EntityPtr   m_Child;              // link this just for the
public:
    Surface( const std::vector< BrushPtr >& assets );

    virtual ~Surface();

protected:
    virtual bool DoInitialize( Renderer* renderer ) throw(std::exception);

    virtual void DoRender() throw(std::exception);

    virtual void DoUpdate( float ticks ) throw(std::exception);

};

#endif /* MESH_H */
