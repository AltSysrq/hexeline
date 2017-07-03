out vec4 dst;

uniform vec4 colour;
uniform sampler2D tex;

in vec2 v_texcoord;

void main() {
  dst = texture2D(tex, v_texcoord);
}
