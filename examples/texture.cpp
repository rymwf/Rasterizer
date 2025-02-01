#include <iostream>
#include <vector>

#include "glm/ext.hpp"
#include "image.h"
#include "pipeline.h"

struct Shader : public ShaderBase<Shader> {
    struct VertexAttribute {
        glm::vec3 pos;
        glm::vec2 uv;
    };

    struct VertOut {
        glm::vec2 uv;  //

        VertOut operator+(VertOut const& other) const {
            VertOut ret = *this;
            ret.uv += other.uv;
            return ret;
        }
        VertOut operator*(float c) const {
            VertOut ret = *this;
            ret.uv *= c;
            return ret;
        }
    };

    struct FragOut {
        glm::vec4 color;
    };

    ImageRGBA8* image;

    glm::mat4 PV;

    PerVertex vert(VertexAttribute const& vertex) {
        PerVertex ret;
        ret.pos = PV * glm::vec4(vertex.pos, 1);
        ret.vertOut.uv = vertex.uv;
        return ret;
    }
    FragOut frag(FragIn const& frag) {
        FragOut ret;
        // ret.color = glm::vec4(frag.vertOut.uv, 0, 1);
        ret.color = image->Sample(frag.vertOut.uv.x, frag.vertOut.uv.y);
        return ret;
    }
};

glm::u8vec4 linear2srgbu8(glm::vec4 const& col) { return glm::pow(col, glm::vec4(0.4545)) * 255.f; }

class Framebuffer : public FramebufferBase {
public:
    Image<glm::u8vec4> color;
    Image<float> depth;

    glm::vec4 clearColor{};
    float clearDepth = 1;

    Framebuffer(int width, int height) : color(width, height), depth(width, height) {}

    void Clear(bool colorBuffer, bool depthBuffer) {
        // if (colorBuffer) color.fill(glm::u8vec4(clearColor * 255.f));
        // enable srgb
        if (colorBuffer) color.fill(glm::u8vec4(glm::pow(clearColor, glm::vec4(0.4545)) * 255.f));
        if (depthBuffer) depth.fill(clearDepth);
    }

    void SetPixel(int i, int j, void const* usrData) override {
        auto fragOut = reinterpret_cast<Shader::FragOut const*>(usrData);
        color(i, j) = linear2srgbu8(glm::clamp(fragOut->color, glm::vec4(0), glm::vec4(1)));
    }
    void SetDepth(int i, int j, float d) override { depth(i, j) = d; }
    float GetDepth(int i, int j) const override { return depth(i, j); }
};

#define PI 3.141592653

struct Mesh {
    PrimitiveTopology topology;
    std::vector<Shader::VertexAttribute> vertices;
};

Mesh testQuad() {
    Mesh mesh;
    mesh.topology = PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;
    mesh.vertices = {
        {{-1, 0, -1}, {0, 0}},  //
        {{1, 0, 1}, {1, 1}},    //
        {{1, 0, -1}, {1, 0}},   //

        {{-1, 0, -1}, {0, 0}},  //
        {{-1, 0, 1}, {0, 1}},   //
        {{1, 0, 1}, {1, 1}},    //
    };
    return mesh;
}

int main() {
    auto mesh = testQuad();

    ImageRGBA8 image;
    image.Load("0.png");

    // config fbo
    Framebuffer fbo(600, 600);
    fbo.clearColor = {0.2, 0.2, 0.2, 1};

    // config pipeline
    Pipeline<Shader, Shader::VertexAttribute> pipeline;
    pipeline.fbo = &fbo;
    pipeline.viewport = {0, 0, 600, 600, -1, 1};
    pipeline.enableDepthTest = true;

    pipeline.shader.image = &image;
    pipeline.shader.PV = glm::perspective(glm::radians(60.f), 1.f, 0.1f, 10.f) *
                         glm::lookAt(glm::vec3(0, 2, 2), glm::vec3(0), glm::vec3(0, 1, 0));

    // clear color buffer
    fbo.Clear(true, true);

    pipeline.Draw(mesh.topology, mesh.vertices);

    // gamma correct
    fbo.color.SavePNG(__FILE__ ".png");
}