#include <iostream>
#include <vector>

#include "image.h"
#include "pipeline.h"

#include "image.h"

struct Shader : public ShaderBase<Shader> {
    struct VertexAttribute {
        glm::vec3 pos;
        glm::vec3 color;
    };

    struct VertOut {
        glm::vec3 col;  //
        VertOut operator+(VertOut const& other) const {
            VertOut ret = *this;
            ret.col += other.col;
            return ret;
        }
        VertOut operator*(float c) const {
            VertOut ret = *this;
            ret.col *= c;
            return ret;
        }
    };

    struct FragOut {
        glm::vec4 color;
    };

    PerVertex vert(VertexAttribute const& vertex) {
        PerVertex ret;
        ret.pos = glm::vec4(vertex.pos, 1);
        ret.vertOut.col = vertex.color;
        return ret;
    }
    FragOut frag(FragIn const& frag) {
        FragOut ret;
        ret.color = glm::vec4(frag.vertOut.col, 1);
        // ret.color = glm::vec4(glm::vec3(frag.fragCoord.z), 1);
        return ret;
    }
};

glm::u8vec4 linear2srgbu8(glm::vec4 const& col) { return glm::pow(col, glm::vec4(0.4545)) * 255.f; }

class Framebuffer : public FramebufferBase {
public:
    ImageRGBA8 color;
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

Mesh testLines0() {
    Mesh mesh;
    mesh.topology = PRIMITIVE_TOPOLOGY_LINE_LIST;
    // mesh.topology = PRIMITIVE_TOPOLOGY_POINT_LIST;
    const int N = 10;
    Shader::VertexAttribute p{};
    for (int i = 0; i < N; ++i) {
        auto a = float(i) / N * 2 * PI;
        auto cosa = std::cos(a);
        auto sina = std::sin(a);
        p.pos = {0, 0, 0.5};
        p.color = {0, 1, 0};
        mesh.vertices.emplace_back(p);
        p.pos = {cosa * 1.0, sina * 1.0, 0.5};
        p.color = {1, 0, 0.0};
        mesh.vertices.emplace_back(p);
    }
    return mesh;
};

Mesh testLines1() {
    Mesh mesh;
    mesh.topology = PRIMITIVE_TOPOLOGY_LINE_LIST;
    mesh.vertices = {
        {{-1, -1, 0}, {1, 0, 0}},
        {{1, -1, 0}, {0, 1, 0}},

        {{1, -1, 0}, {0, 1, 0}},
        {{1, 1, 0}, {1, 0, 0}},

        {{1, 1, 0}, {1, 0, 0}},
        {{-1, 1, 0}, {0, 1, 0}},

        {{-1, 1, 0}, {0, 1, 0}},
        {{-1, -1, 0}, {1, 0, 0}},
    };
    return mesh;
}
Mesh testLines2() {
    Mesh mesh;
    mesh.topology = PRIMITIVE_TOPOLOGY_LINE_LIST;
    mesh.vertices = {
        // {{-2, 0, 0}, {1, 0, 0}},//
        // {{2, 0, 0}, {0, 1, 0}},//

        // {{0, -2, 0}, {1, 0, 0}},//
        // {{0, 2, 0}, {0, 1, 0}},//

        {{-1, 0.5,-2}, {1, 0, 0}},//
        {{1, 0.5,2}, {0, 1, 0}},//
    };
    return mesh;
}

Mesh testTriangle0() {
    Mesh mesh;
    mesh.topology = PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;
    mesh.vertices = {
        {{-1.5, -0.5, 0.1}, {1, 0, 0}},
        {{1.5, -0.8, 0.5}, {0, 1, 0}},
        {{0, 1.3, 2.0}, {0, 0, 1}},
    };
    return mesh;
}
Mesh testTriangle1() {
    Mesh mesh;
    mesh.topology = PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;
    mesh.vertices = {
        {{-0.5, -0.5, 0.1}, {1, 0, 0}},
        {{2, -1, 0.8}, {0, 1, 0}},
        {{0, 1.5, 0.2}, {0, 0, 1}},

        {{-2, -1, 0.8}, {1, 0, 0}},
        {{0.5, 0.1, 0.1}, {0, 1, 0}},
        {{-1, 1.5, 0.5}, {0, 0, 1}},
    };
    return mesh;
}
Mesh testTriangle2() {
    Mesh mesh;
    mesh.topology = PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;
    mesh.vertices = {
        {{0,0,0},{1,0,0}},
        {{1,0,0},{0,1,0}},
        {{0,1,0},{0,0,1}},
    };
    return mesh;
}

int main() {
    // auto mesh = testLines0();
    // auto mesh = testLines1();
    // auto mesh = testLines2();
    // auto mesh = testTriangle0();
    auto mesh = testTriangle1();
    // auto mesh = testTriangle2();

    // config fbo
    Framebuffer fbo(600, 600);
    fbo.clearColor = {0.2, 0.2, 0.2, 1};

    // config pipeline
    Pipeline<Shader, Shader::VertexAttribute> pipeline;
    pipeline.fbo = &fbo;
    pipeline.viewport = {0, 0, 600, 600, -1, 1};
    pipeline.enableDepthTest = true;

    // clear color buffer
    fbo.Clear(true, true);

    pipeline.Draw(mesh.topology, mesh.vertices);

    // gamma correct
    fbo.color.SavePNG(__FILE__ ".png");
}