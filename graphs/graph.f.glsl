#version 430

in vec3 LightIntensity;

layout( location = 0 ) out vec4 FragColor;

uniform int switch_color;
uniform vec4 color;
uniform vec4 colortest;
uniform vec4 colortext;
varying vec2 texpos;
uniform sampler2D tex;
varying vec4 graph_coord;

void main()
{
	float factor;
	if (gl_FrontFacing)
		factor = 1.0;
	else
		factor = 0.5;
	if ( switch_color == 1 )
	{
		FragColor = color;
	}
	else
	{
		if ( switch_color == 2 )
			FragColor = vec4(1, 1, 1, texture2D(tex, texpos).a) * colortext;
		else
		{
			if ( switch_color == 3 )
				FragColor = colortest;
			else
				FragColor = ( ( /* graph_coord * 0.0 */color * 0.0 + graph_coord.z ) / 1.0 + 0.0 ) * factor;
		}
	}
}
