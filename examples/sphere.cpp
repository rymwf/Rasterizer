#include <iostream>
#include <vector>

#include "glm/ext.hpp"
#include "image.h"
#include "pipeline.h"

struct Shader : public ShaderBase<Shader> {
    struct VertexAttribute {
        glm::vec3 pos;
        glm::vec3 normal;
    };

    struct VertOut {
        glm::vec3 normal;  //

        VertOut operator+(VertOut const& other) const {
            VertOut ret = *this;
            ret.normal += other.normal;
            return ret;
        }
        VertOut operator*(float c) const {
            VertOut ret = *this;
            ret.normal *= c;
            return ret;
        }
    };

    struct FragOut {
        glm::vec4 color;
    };

    glm::mat4 P;
    glm::mat4 V;

    PerVertex vert(VertexAttribute const& vertex) {
        PerVertex ret;
        ret.pos = P * V * glm::vec4(vertex.pos, 1);
        ret.vertOut.normal = glm::vec3(V * glm::vec4(vertex.normal, 0.f));
        return ret;
    }
    FragOut frag(FragIn const& frag) {
        FragOut ret;

        const auto L = glm::normalize(glm::vec3(1, 1, 1));

        auto NdotL = glm::dot(frag.vertOut.normal, L);
        auto col = glm ::vec3(NdotL + 0.2);
        ret.color = glm::vec4(col, 1);
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
    std::vector<unsigned int> indices;
};

Mesh CreateSphereMesh(float radius, int longitudeSegCount, int latitudeSegCount) {
    Mesh mesh;
    mesh.topology = PRIMITIVE_TOPOLOGY_TRIANGLE_STRIP;
    // mesh.topology = PRIMITIVE_TOPOLOGY_LINE_STRIP;

    float y, ry;
    float theta, phi, u, v;
    Shader::VertexAttribute vertex;
    for (int j = 0; j <= longitudeSegCount; ++j) {
        v = float(j) / longitudeSegCount;
        theta = v * 3.14159265359;
        y = -std::cos(theta);
        ry = std::sin(theta);
        for (int i = 0; i < latitudeSegCount; ++i) {
            u = float(i) / latitudeSegCount;
            phi = u * 6.28318530718;
            vertex.normal = {ry * std::cos(phi), y, ry * std::sin(phi)};
            vertex.pos = vertex.normal * radius;
            mesh.vertices.emplace_back(vertex);
        }
    }
    mesh.vertices.shrink_to_fit();
    mesh.indices.reserve((2 * (latitudeSegCount + 1) + 1) * longitudeSegCount);

    for (int j = 0; j < longitudeSegCount; ++j) {
        for (int i = 0; i < latitudeSegCount; ++i) {
            mesh.indices.emplace_back(j * latitudeSegCount + i);
            mesh.indices.emplace_back((j + 1) * latitudeSegCount + i);
        }
        mesh.indices.emplace_back(j * latitudeSegCount);
        mesh.indices.emplace_back((j + 1) * latitudeSegCount);
        mesh.indices.emplace_back(~0u);
    }
    return mesh;
}

int main() {
    auto mesh = CreateSphereMesh(1, 32, 32);

    // config fbo
    Framebuffer fbo(600, 600);
    fbo.clearColor = {0.2, 0.2, 0.2, 1};

    // config pipeline
    Pipeline<Shader, Shader::VertexAttribute> pipeline;
    pipeline.fbo = &fbo;
    pipeline.viewport = {0, 0, 600, 600, -1, 1};
    pipeline.enableDepthTest = true;

    pipeline.shader.P = glm::perspective(glm::radians(60.f), 1.f, 0.1f, 10.f);
    pipeline.shader.V = glm::lookAt(glm::vec3(0, 0, 3), glm::vec3(0), glm::vec3(0, 1, 0));

    // clear color buffer
    fbo.Clear(true, true);

    pipeline.DrawIndices(mesh.topology, mesh.vertices, mesh.indices);

    // gamma correct
    fbo.color.SavePNG(__FILE__ ".png");
}