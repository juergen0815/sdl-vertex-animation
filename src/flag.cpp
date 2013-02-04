/*
 * cube.cpp
 *
 *  Created on: Dec 4, 2011
 *      Author: Jurgen
 */


#include "flag.h"

#include <cmath>

const int columns = 44;
const int rows   = 44;

struct PoolDeleter
{
    void operator()(void const *p) const
    {
        // called when shared_ptr goes out of scope - that is last container releases link to pool
        EntityPool::DestroyPool( (MemoryPool*)p );
    }
};

Flag::Flag( const std::vector< BrushPtr >& assets )
    : m_VboID(-1)
    , m_Assets( assets )
    , m_MemoryPool( EntityPool::CreatePool<Vector>( 0 ), PoolDeleter() )
    , m_VertexBuffer( m_MemoryPool )    // use the same memory pool for vertex and texture coords
    , m_TexCoordBuffer( m_MemoryPool )
{
    // we might just want to create this in DoInitialize - and throw away the data we don't need locally

    // allocate memory buffers for vertex and texture coord
    m_VertexBuffer.resize( columns*rows );
    m_TexCoordBuffer.resize( columns*rows );

    // generate vertex array
    float vx;
    float vy;
    float vz;
    const float step = 0.2f;
    const float amp  = 2.0f;
    const float loop = 8.0f; // num of sin loops (or waves)
    auto vertex = m_VertexBuffer.begin();
    for ( float y = 0; y < rows; ++y )
    {
        for ( float x = 0; x < columns; ++x, ++vertex )
        {
            vertex[ Vector::X ] = x * step - 4.4f;
            vertex[ Vector::Y ] = y * step - 4.4f;
            // maybe I should shift this for each row, huh, norm x to "length" of column (0.0 - 1.0)
            vertex[ Vector::Z ] = std::sin( (x / columns) * ( M_PI / 180.0f * loop) ) * amp; // make z a big "wavy"
        }
    }
    // generate index array; we got rows * columns * 2 tris
    m_IndexArray.resize( rows * columns * 3 * 2 ); // 3 vertices per tri, 2 tri per quad = 6 entries per iteration
    int looper(0);
    for ( ; looper < rows*columns; ++looper ) {
        // map 4 quad vertices to 6 tris
        // top tri 0 - 1 - 2
        m_IndexArray[ looper + 0 ] = 0+looper;
        m_IndexArray[ looper + 1 ] = 1+looper;
        m_IndexArray[ looper + 2 ] = 3+looper;
        // bottom tri 3- 1- 2
        m_IndexArray[ looper + 3 ] = 3+looper;
        m_IndexArray[ looper + 4 ] = 1+looper;
        m_IndexArray[ looper + 5 ] = 2+looper;
    }
}

Flag::~Flag()
{

}

bool Flag::DoInitialize( Renderer* renderer ) throw(std::exception)
{
    bool r(false);

    bool hasVBO  = glewGetExtension("GL_ARB_vertex_buffer_object");
    ASSERT( hasVBO, "VBOs not supported!" );
    if ( hasVBO ) {
        glGenBuffers(1, &m_VboID);
        glBindBuffer(GL_ARRAY_BUFFER, m_VboID);
        glBufferData(GL_ARRAY_BUFFER, sizeof(Vector)*m_VertexBuffer.size(), 0, GL_STATIC_DRAW_ARB);
        std::size_t offset(0);
        // copy vertices starting from 0 offest
        float *vertices = (float*)&m_VertexBuffer[0];
        glBufferSubData(GL_ARRAY_BUFFER, offset, sizeof(Vector)*m_VertexBuffer.size(), vertices);
        offset += sizeof(Vector)*m_VertexBuffer.size();
        glBufferSubData(GL_ARRAY_BUFFER, offset, sizeof(int)*m_IndexArray.size(), &m_IndexArray[0] );
#if 0
        offset += sizeof(vertices);
        glBufferSubData(GL_ARRAY_BUFFER, offset, sizeof(normals), normals);                // copy normals after vertices
        offset += sizeof(normals);
        glBufferSubData(GL_ARRAY_BUFFER, offset, sizeof(colors), colors);  // copy colors after normals
        if ( m_Texture ) {
            offset += sizeof(texCoords);
            glBufferSubData(GL_ARRAY_BUFFER, offset, sizeof(texCoords), texCoords);  // copy colors after normals
        }
#endif
    }
    r = true;

    return r;
}

void Flag::DoRender() throw(std::exception)
{
#if 0
    // enable vertex arrays
    int normalArrayEnabled;
    glGetIntegerv( GL_NORMAL_ARRAY, &normalArrayEnabled );
    if ( !normalArrayEnabled )  {
        glEnableClientState(GL_NORMAL_ARRAY);
    }
    int colorArrayEnabled;
    glGetIntegerv( GL_COLOR_ARRAY, &colorArrayEnabled );
    if ( !colorArrayEnabled && (GetRenderState()->GetFlags() & BLEND_COLOR_F) ) {
        glEnableClientState(GL_COLOR_ARRAY);
    }
#endif
    int vertexArrayEnabled;
    glGetIntegerv( GL_VERTEX_ARRAY, &vertexArrayEnabled );
    if (!vertexArrayEnabled) {
        glEnableClientState(GL_VERTEX_ARRAY);
    }
    int m_TexCoordArrayEnabled;
    if ( m_Texture && !m_TexCoordArrayEnabled ) {
        glGetIntegerv( GL_TEXTURE_COORD_ARRAY, &m_TexCoordArrayEnabled );
        if (!m_TexCoordArrayEnabled) {
            glEnableClientState(GL_TEXTURE_COORD_ARRAY);
        }
    }

    int blend_enabled;
    glGetIntegerv(GL_BLEND, &blend_enabled);

#ifdef _USE_VBO_
    glBindBuffer(GL_ARRAY_BUFFER, m_VboID);
    // before draw, specify vertex and index arrays with their offsets
    std::size_t offset(0);
    glVertexPointer(3, GL_FLOAT, 0, (void*)offset); offset += sizeof(vertices);
    glNormalPointer(   GL_FLOAT, 0, (void*)offset); offset += sizeof(normals);
    glColorPointer (4, GL_FLOAT, 0, (void*)offset);
    if ( m_Texture) {
        m_Texture->Enable();
        offset += sizeof(texCoords);
        glTexCoordPointer(2, GL_FLOAT, 0, (void*)offset );
        if (blend_enabled) {
            glDisable( GL_BLEND );
        }
    }
    // use draw indices
    glDrawArrays(GL_TRIANGLES, 0, 36);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
#else
        // before draw, specify vertex arrays
        float *vertices = (float*)&m_VertexBuffer[0];
        glVertexPointer(4, GL_FLOAT, 0, vertices);
#if 0
        glNormalPointer(  GL_FLOAT, 0, normals);
        glColorPointer(4, GL_FLOAT, 0, colors);
#endif
        if ( m_Texture) {
            m_Texture->Enable();
            float *texCoords = (float*)&m_TexCoordBuffer[0];
            glTexCoordPointer( 4, GL_FLOAT, 0, texCoords );
        }
        int *indices = &m_IndexArray[0];
        glDrawElements(GL_TRIANGLES, m_IndexArray.size(), GL_UNSIGNED_INT, indices );
#endif

    if (!vertexArrayEnabled)  {
        glDisableClientState(GL_VERTEX_ARRAY);  // disable vertex arrays
    }
#if 0
    if (!colorArrayEnabled)   {
        glDisableClientState(GL_COLOR_ARRAY);
    }
    if (!normalArrayEnabled) {
        glDisableClientState(GL_NORMAL_ARRAY);
    }
#endif
    if ( m_Texture ) {
        m_Texture->Disable();
        if (!m_TexCoordArrayEnabled) {
            glDisableClientState(GL_TEXTURE_COORD_ARRAY);
        }
        if (blend_enabled) {
            glEnable( GL_BLEND );
        }
    }


}

void Flag::DoUpdate( float ticks ) throw(std::exception)
{
    ASSERT( m_VertexBuffer <= rows*columns, "Vertex array to small!" );

    // TODO: Do vertex animation
    auto vertex = m_VertexBuffer.begin();
    for ( float y = 0; y < rows; ++y )
    {
        for ( float x = 0; x < columns; ++x, ++vertex )
        {
//            vx = x * step;
//            vy = y * step;
//            vz = std::sin( x / colums * M_PI / 360.0f ) * amp; // make z a big "wavy"
//            vertex[ Vector::X ] = vx;
//            vertex[ Vector::Y ] = vy;
//            vertex[ Vector::Z ] = vz;
        }
    }
}

#if 0
#include "config.h"
#include "wave_mesh.h"

#include <render/tex_surface.h>

namespace wave {

	static GLfloat points[45][45][3];    // the array for the points on the grid of our "wave"

	class VertexTarget : public pei::render::VertexArrayRenderTarget
	{
	public:
		void AllocateVertexBuffer( unsigned int num )
		{
			m_Vertices.resize( num );
		}
		GLfloat* GetVertexBuffer() { return &m_Vertices.at(0); }

		GLfloat* GetTexCoordBuffer() { return &m_TextureCoords.at(0); }
	};

	WaveMesh::WaveMesh()
		: amp(4)
		, xrot(0)
    	, yrot()
    	, zrot()
		, m_TimeEllapsed(0)
		, m_Speed(0.75)
	{
		VertexTarget *rt = new VertexTarget;
		SetRenderTarget( pei::RenderTargetInterfacePtr(rt) );
	    float float_x, float_y; // loop counters.

	    for(float_x = 0.0f; float_x < 9.0f; float_x +=  0.2f )	{
			for(float_y = 0.0f; float_y < 9.0f; float_y += 0.2f)		{
				points[ (int) (float_x*5) ][ (int) (float_y*5) ][0] = float_x - 4.4f;
				points[ (int) (float_x*5) ][ (int) (float_y*5) ][1] = float_y - 4.4f;
				points[ (int) (float_x*5) ][ (int) (float_y*5) ][2] = (float) (sin( ( (float_x*5*8)/360 ) * 3.14159 * amp ));
			}
	    }

		// indices and texCoords remain static, only vertices will be updated
	    std::vector<unsigned int> indexArray;
	    int indexOffset(0);
	    int x, y;
		unsigned int numVertices(0);
	    float float_xb, float_yb;
	    for (x=0; x<44; x++) {
			for (y=0; y<44; y++) {
				float_x  = (float) (x)/44;
				float_y  = (float) (y)/44;
				float_xb = (float) (x+1)/44;
				float_yb = (float) (y+1)/44;

				texCoords.push_back(float_x );  texCoords.push_back( float_y);
				texCoords.push_back(float_x );  texCoords.push_back( float_yb);
				texCoords.push_back(float_xb ); texCoords.push_back( float_yb);
				texCoords.push_back(float_xb ); texCoords.push_back( float_y);

				// must restart poly for proper texture mapping...optimize this
				indexArray.push_back(0+numVertices);
				indexArray.push_back(1+numVertices);
				indexArray.push_back(3+numVertices);
				indexArray.push_back(3+numVertices);
				indexArray.push_back(1+numVertices);
				indexArray.push_back(2+numVertices);
				// use TRIANGLES, cannot use TRI_STRIPS - triangle getting merged into 1 draw call, strips don't
				rt->AddPolygon( GL_TRIANGLES, 6, indexOffset );
				indexOffset += 6;
				numVertices += 4;
			}
	    }
		rt->SetTexCoords(texCoords.size(), &texCoords.at(0) );
		rt->SetIndices( indexArray.size(), &indexArray.at(0));
		rt->AllocateVertexBuffer( numVertices*3 );
		SetColor( pei::Color( 1.0,1.0,1.0 ));
		RequestStateChanged( true );
	}

	WaveMesh::~WaveMesh()
	{
	}

	pei::SurfacePtr WaveMesh::SetSurface(pei::SurfacePtr s, int n)
	{
		pei::Drawable::SetSurface( s, n );
#if 0
		// read back texture
		pei::Texture2DSurfacePtr tex = boost::dynamic_pointer_cast< pei::Texture2DSurface >( GetSurface(0) );
		if ( tex ) {
			// We must re calc texture coords for each "frame" we attach to the cube - for non pow2 surfaces!
//			int numCoords = texCoords.size();
			VertexTarget *rt = (VertexTarget*)(GetRenderTarget().get());
			GLfloat *texCoords = rt->GetTexCoordBuffer();
			const pei::TexSurface::TexCoords &triCoords = tex->GetCoords();

			int n(0);
			float x,y;
		    float tx0, tx1, ty0, ty1;
		    float dx = (triCoords.MaxX-triCoords.MinX)/44;
		    float dy = (triCoords.MaxY-triCoords.MinY)/44;

		    for (x = triCoords.MinX; x < triCoords.MaxX; x += dx ) {
				for (y = triCoords.MinY; y < triCoords.MaxY; y += dy) {
					tx0 = x;
					ty0 = y;
					tx1 = x + dx;
					ty1 = y + dy;

					texCoords[ n++ ] = tx0; texCoords[ n++ ] = ty0;
					texCoords[ n++ ] = tx0; texCoords[ n++ ] = ty1;
					texCoords[ n++ ] = tx1; texCoords[ n++ ] = ty1;
					texCoords[ n++ ] = tx1; texCoords[ n++ ] = ty0;
				}
		    }
		}
#endif
		return s;
	}

	void WaveMesh::OnUpdateAnimation( const pei::RenderProfilePtr& profile, double time )
	{
		m_TimeEllapsed += time;
	    if ( m_TimeEllapsed*m_Speed > 16.67 ) {
	    	m_TimeEllapsed = 0;
			int x, y;
			for (y = 0; y <45; y++) {
				points[44][y][2] = points[0][y][2];
			}

			for( x = 0; x < 44; x++ ) {
				for( y = 0; y < 45; y++) {
					points[x][y][2] = points[x+1][y][2];
				}
			}
			// no point in updating if vertices don't change
			RequestStateChanged( true );

//		    xrot +=0.3f;
//		    yrot +=0.2f;
//		    zrot +=0.4f;
		    SetRotation( GetXRotation(), GetYRotation() + 0.2f , GetZRotation() );
	    }
	}

	void WaveMesh::OnUpdateState( const pei::RenderProfilePtr& profile )
	{
		if ( MustRefreshState( ) )
		{
			VertexTarget *rt = (VertexTarget*)(GetRenderTarget().get());
			GLfloat *vertexArray = rt->GetVertexBuffer();
			int x, y;
			int idx(0);
			for (x=0; x<44; x++) {
				for (y=0; y<44; y++) {
				    // 0
					vertexArray[idx++] = ( points[x][y][0] );
					vertexArray[idx++] = ( points[x][y][1] );
					vertexArray[idx++] = ( points[x][y][2] );

					// 1
					vertexArray[idx++] = ( points[x][y+1][0] );
					vertexArray[idx++] = ( points[x][y+1][1] );
					vertexArray[idx++] = ( points[x][y+1][2] );

					// 2
					vertexArray[idx++] = ( points[x+1][y+1][0] );
					vertexArray[idx++] = ( points[x+1][y+1][1] );
					vertexArray[idx++] = ( points[x+1][y+1][2] );

					// 3
					vertexArray[idx++] = ( points[x+1][y][0] );
					vertexArray[idx++] = ( points[x+1][y][1] );
					vertexArray[idx++] = ( points[x+1][y][2] );
				}
			}
		}
	}

	void WaveMesh::OnDraw( const pei::RenderProfilePtr& profile, const pei::SurfacePtr& buffer, const pei::RenderParam& param  )
	{
		pei::Drawable::OnDraw( profile, buffer, param );
	}

}
#endif
