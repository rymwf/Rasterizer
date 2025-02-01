#include <iostream>
#include <vector>
#include <filesystem>

#include "glm/ext.hpp"
#include "image.h"
#include "pipeline.h"

#include "image.h"

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

struct Shader : public ShaderBase<Shader> {
    struct VertexAttribute {
        glm::vec3 pos;
        glm::vec3 normal;
        glm::vec2 uv;
    };

    struct VertOut {
        glm::vec3 normal;  //
        glm::vec2 uv;  //

        VertOut operator+(VertOut const& other) const {
            VertOut ret = *this;
            ret.normal+= other.normal;
            ret.uv += other.uv;
            return ret;
        }
        VertOut operator*(float c) const {
            VertOut ret = *this;
            ret.normal *= c;
            ret.uv *= c;
            return ret;
        }
    };

    struct FragOut {
        glm::vec4 color;
    };

    glm::mat4 P;
    glm::mat4 V;
    ImageRGBA8* image;
    ImageRGBA8* imageEmission;

    PerVertex vert(VertexAttribute const& vertex) {
        PerVertex ret;
        ret.pos = P * V * glm::vec4(vertex.pos, 1);
        ret.vertOut.normal = glm::vec3(V * glm::vec4(vertex.normal, 0.f));
        ret.vertOut.uv = vertex.uv;
        return ret;
    }
    FragOut frag(FragIn const& frag) {
        FragOut ret;

        const auto L = glm::normalize(glm::vec3(1, 1, 1));

        auto NdotL = glm::dot(frag.vertOut.normal, L);
        auto col = (NdotL + 0.2f) * image->Sample(frag.vertOut.uv.x, frag.vertOut.uv.y);
        col += imageEmission->Sample(frag.vertOut.uv.x, frag.vertOut.uv.y);
        ret.color = col;
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

    ImageRGBA8 imageDiff;
    ImageRGBA8 imageEmission;
};

Mesh LoadObj(std::string const& inputfile) {
    Mesh mesh;
    mesh.topology = PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;

    tinyobj::ObjReaderConfig reader_config;
    std::filesystem::path path(inputfile);
    reader_config.mtl_search_path = path.parent_path().string();

    tinyobj::ObjReader reader;

    if (!reader.ParseFromFile(inputfile, reader_config)) {
        if (!reader.Error().empty()) {
            std::cerr << "TinyObjReader: " << reader.Error();
        }
        exit(1);
    }

    if (!reader.Warning().empty()) {
        std::cout << "TinyObjReader: " << reader.Warning();
    }

    auto& attrib = reader.GetAttrib();
    auto& shapes = reader.GetShapes();
    auto& materials = reader.GetMaterials();

    // mesh.diffTexName= +materials[0].diffuse_texname;
    mesh.imageDiff.Load((reader_config.mtl_search_path + "/" + materials[0].diffuse_texname).c_str());
    mesh.imageEmission.Load((reader_config.mtl_search_path + "/" + materials[0].emissive_texname).c_str());

    // Loop over shapes

    Shader::VertexAttribute vertex;

    for (size_t s = 0; s < shapes.size(); s++) {
        // Loop over faces(polygon)
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
            size_t fv = size_t(shapes[s].mesh.num_face_vertices[f]);

            // Loop over vertices in the face.
            for (size_t v = 0; v < fv; v++) {
                // access to vertex
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
                tinyobj::real_t vx = attrib.vertices[3 * size_t(idx.vertex_index) + 0];
                tinyobj::real_t vy = attrib.vertices[3 * size_t(idx.vertex_index) + 1];
                tinyobj::real_t vz = attrib.vertices[3 * size_t(idx.vertex_index) + 2];

                vertex.pos = {vx, vy, vz};

                // Check if `normal_index` is zero or positive. negative = no normal data
                if (idx.normal_index >= 0) {
                    tinyobj::real_t nx = attrib.normals[3 * size_t(idx.normal_index) + 0];
                    tinyobj::real_t ny = attrib.normals[3 * size_t(idx.normal_index) + 1];
                    tinyobj::real_t nz = attrib.normals[3 * size_t(idx.normal_index) + 2];
                    vertex.normal = {nx, ny, nz};
                }

                // Check if `texcoord_index` is zero or positive. negative = no texcoord data
                if (idx.texcoord_index >= 0) {
                    tinyobj::real_t tx = attrib.texcoords[2 * size_t(idx.texcoord_index) + 0];
                    tinyobj::real_t ty = attrib.texcoords[2 * size_t(idx.texcoord_index) + 1];
                    vertex.uv = {tx, ty};
                }

                // Optional: vertex colors
                // tinyobj::real_t red   = attrib.colors[3*size_t(idx.vertex_index)+0];
                // tinyobj::real_t green = attrib.colors[3*size_t(idx.vertex_index)+1];
                // tinyobj::real_t blue  = attrib.colors[3*size_t(idx.vertex_index)+2];
                mesh.vertices.emplace_back(vertex);
            }
            index_offset += fv;

            // per-face material
            shapes[s].mesh.material_ids[f];
        }
    }

    return mesh;
}

int main() {
    auto mesh = LoadObj("DamagedHelmet/DamagedHelmet.obj");

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
    pipeline.shader.image = &mesh.imageDiff;
    pipeline.shader.imageEmission = &mesh.imageEmission;

    // clear color buffer
    fbo.Clear(true, true);

    pipeline.Draw(mesh.topology, mesh.vertices);

    // gamma correct
    fbo.color.SavePNG(__FILE__ ".png");
}