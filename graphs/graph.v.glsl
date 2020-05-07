#version 430

layout (location = 0) in vec3 VertexPosition;
layout (location = 1) in vec3 VertexNormal;

out vec3 LightIntensity;

attribute vec2 coord2d;
attribute vec4 coordtext;
varying vec4 graph_coord;
uniform int switch_transform;
uniform mat4 transform;
uniform mat4 texture_transform;
uniform mat4 vertex_transform;
uniform mat4 vertex_transform_textX;
uniform mat4 vertex_transform2; // 3
uniform mat4 vertex_transform90; // 4
uniform mat4 vertex_transform180_90; // 5
uniform mat4 vertex_transform180_90_vert; // 6
uniform sampler2D mytexture;
varying vec2 texpos;

void main()
{
	graph_coord = texture_transform * vec4(coord2d, 0, 1);
	graph_coord.z = (texture2D(mytexture, graph_coord.xy / 2.0 + 0.5).r);
	if ( switch_transform == 1 )
	{
		gl_Position = vertex_transform * vec4(coord2d, 0.0, 1);
	}
	else
	{
		if ( switch_transform == 2 )
		{
			gl_Position = vertex_transform * vec4(coordtext.xy, 0, 1);
			texpos = coordtext.zw;
		}
		else
		{
			if ( switch_transform == 3 )
				gl_Position = vertex_transform2 * vec4(coord2d, 0.0, 1);
			else
			{
				if ( switch_transform == 4 )
					gl_Position = vertex_transform90 * vec4(coord2d, 0.0, 1);
				else
				{
					if ( switch_transform == 5 )
						gl_Position = vertex_transform180_90 * vec4(coord2d, 0.0, 1);
					else
					{
						if ( switch_transform == 6 )
							gl_Position = vertex_transform180_90_vert * vec4(coord2d, 0.0, 1);
						else
						{
							if ( switch_transform == 7 )
							{
								gl_Position = vertex_transform_textX * vec4(coordtext.xy, 0, 1);
								texpos = coordtext.zw;
							}
							else
								gl_Position = vertex_transform * vec4(coord2d, graph_coord.z, 1);
						}
					}
				}
			}
		}
	}

	LightIntensity = vec3(0.7, 0.7, 0.7);
}
