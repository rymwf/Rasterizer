#pragma once
#include <concepts>
#include <optional>
#include <vector>

// template <typename T>
// concept con_shader_t = requires {
//     typename T::VertexAttribute;
//     typename T::VertOut;
//     typename T::FragIn;
//     typename T::FragOut;
//     // { typename T::vert(T::VertexAttribute) } -> std::same_as<typename T::VertOut>;
//     // { typename T::frag(T::FragIn) } -> std::same_as<typename T::FragOut>;
//     //  typename T::frag(T::FragIn const&);// } -> std::same_as<typename T::FragOut>;
// };

struct FragCoord {
    int x;
    int y;
    float z;
    float w;

    FragCoord operator*(float c) const { return FragCoord(x * c, y * c, z * c, w * c); }
    FragCoord operator+(FragCoord const& other) const {
        return FragCoord(x + other.x, y + other.y, z + other.z, w + other.w);
    }
};

template <typename T>
struct ShaderBase {
    struct PerVertex {
        glm::vec4 pos;  // clip space
        typename T::VertOut vertOut;
    };
    struct FragIn {
        FragCoord fragCoord;
        typename T::VertOut vertOut;
    };
    float fragDepth;
};

struct FramebufferBase {
    // int width;
    // int height;
    FramebufferBase() = default;
    // FramebufferBase(int width, int height) : width(width), height(height) {}
    virtual ~FramebufferBase() {}

    virtual void SetPixel(int i, int j, void const* usrData) = 0;

    virtual void SetDepth(int i, int j, float d) {}
    virtual float GetDepth(int i, int j) const { return 0; }
};

//
enum PrimitiveTopology {
    PRIMITIVE_TOPOLOGY_POINT_LIST,
    PRIMITIVE_TOPOLOGY_LINE_LIST,
    PRIMITIVE_TOPOLOGY_LINE_LOOP,
    PRIMITIVE_TOPOLOGY_LINE_STRIP,
    PRIMITIVE_TOPOLOGY_TRIANGLE_LIST,
    PRIMITIVE_TOPOLOGY_TRIANGLE_STRIP,
    PRIMITIVE_TOPOLOGY_TRIANGLE_FAN,
    PRIMITIVE_TOPOLOGY_PATCHES

};

// template <typename T>
// struct Primitive {
//     using vertex_type = T;
//     PrimitiveTopology topology;
//     std::vector<T> vertices;
// };

struct Viewport {
    float x;
    float y;
    float width;
    float height;
    float n{-1};
    float f{1};
};

struct Scissor {
    float x{};
    float y{};
    float width{1};
    float height{1};
};

template <typename T>
T lerp(T const& a, T const& b, float c) {
    return a * (1 - c) + b * c;
}

template <typename Shader_T, typename VertexAttributes_T>
struct Pipeline {
    using shader_t = Shader_T;
    using vertex_attributes_t = VertexAttributes_T;
    using frag_in_t = typename shader_t::FragIn;
    using frag_out_t = typename shader_t::FragOut;
    using per_vertex_t = typename shader_t::PerVertex;

    shader_t shader;
    FramebufferBase* fbo;
    Viewport viewport;
    Scissor scissor{};

    bool isPerspective = false;
    bool enableDepthTest = false;

    bool enableCullFace = true;  // cull back cw 

    void Shading(frag_in_t const& fragIn) {
        shader.fragDepth = fragIn.fragCoord.z;
        if (enableDepthTest) {
            if (shader.fragDepth < fbo->GetDepth(fragIn.fragCoord.x, fragIn.fragCoord.y)) {
                fbo->SetDepth(fragIn.fragCoord.x, fragIn.fragCoord.y, shader.fragDepth);
                auto fragOut = shader.frag(fragIn);
                fbo->SetPixel(fragIn.fragCoord.x, fragIn.fragCoord.y, &fragOut);
            }
        } else {
            auto fragOut = shader.frag(fragIn);
            fbo->SetPixel(fragIn.fragCoord.x, fragIn.fragCoord.y, &fragOut);
        }
    }

    frag_in_t interpolate(frag_in_t const& a, frag_in_t const& b, float c) const {
        frag_in_t ret;
        if (isPerspective) {
            auto y0 = a.vertOut * b.fragCoord.z;
            auto y1 = b.vertOut * a.fragCoord.z;
            ret.fragCoord.z = a.fragCoord.z * b.fragCoord.z / lerp(b.fragCoord.z, a.fragCoord.z, c);
            ret.vertOut = lerp(y0, y1, c) * ret.fragCoord.z;
        } else {
            ret.fragCoord = lerp(a.fragCoord, b.fragCoord, c);
            ret.vertOut = lerp(a.vertOut, b.vertOut, c);
        }
        return ret;
    }

    static glm::vec4 Clip2NDC(glm::vec4 const& p) {
        auto a = 1. / p.w;
        return glm::vec4(p.x * a, p.y * a, p.z * a, p.w);
    }

    // depth range [0, 1]
    glm::mat4 NDC2ViewportMatrix() const {
        return glm::mat4(viewport.width / 2, 0, 0, 0,   // col 0
                         0, viewport.height / 2, 0, 0,  //
                         0, 0, 0.5, 0.0,                //
                         viewport.width / 2 + viewport.x, viewport.height / 2 + viewport.y, 0.5, 1);
    }
    FragCoord CalcFragCoord(glm::vec4 const& clipPos) const {
        auto ndcPos = Clip2NDC(clipPos);
        auto viewportPos = NDC2ViewportMatrix() * glm::vec4(ndcPos.x, ndcPos.y, ndcPos.z, 1);
        // pixel center convention is integer
        viewportPos += glm::vec4(-0.5, -0.5, 0, 0);
        viewportPos.w = 1. / ndcPos.w;  // move to ndc
        return FragCoord{(int)viewportPos.x, (int)viewportPos.y, viewportPos.z, viewportPos.w};
    }

    void RasterizePoint(frag_in_t const& a) { Shading(a); }
    void RasterizeLineHigh(frag_in_t const& a, frag_in_t const& b) {
        int dx = b.fragCoord.x - a.fragCoord.x;
        int dy = b.fragCoord.y - a.fragCoord.y;

        int sx = 1;
        int sy = 1;
        if (dx < 0) {
            sx = -1;
            dx = -dx;
        }
        if (dy < 0) sy = -1;

        int d = 2 * dx - dy;
        for (int x = a.fragCoord.x, y = a.fragCoord.y; y != b.fragCoord.y; y += sy) {
            auto c = std::abs(float(y - a.fragCoord.y) / dy);
            auto fragIn = interpolate(a, b, c);
            fragIn.fragCoord.x = x;
            fragIn.fragCoord.y = y;

            Shading(fragIn);

            if (d > 0) {
                x += sx;
                d -= 2 * dy;
            }
            d += 2 * dx;
        }
    }
    void RasterizeLineLow(frag_in_t const& a, frag_in_t const& b) {
        //
        int dx = b.fragCoord.x - a.fragCoord.x;
        int dy = b.fragCoord.y - a.fragCoord.y;

        int sx = 1;
        int sy = 1;
        if (dy < 0) {
            sy = -1;
            dy = -dy;
        }
        if (dx < 0) sx = -1;

        int d = 2 * dy - dx;

        for (int x = a.fragCoord.x, y = a.fragCoord.y; x != b.fragCoord.x; x += sx) {
            auto c = std::abs(float(x - a.fragCoord.x) / dx);
            auto fragIn = interpolate(a, b, c);
            fragIn.fragCoord.x = x;
            fragIn.fragCoord.y = y;

            Shading(fragIn);

            if (d > 0) {
                y += sy;
                d -= 2 * dx;
            }
            d += 2 * dy;
        }
    }

    // Function to check if a point is inside the triangle using edge functions
    static bool IsInsideTriangle(glm::ivec2 const& p, glm::ivec2 const& v0, glm::ivec2 const& v1,
                                 glm::ivec2 const& v2) {
        auto cross = [](glm::ivec2 const& e0, glm::ivec2 const& e1) { return e0.x * e1.y - e0.y * e1.x; };

        auto pv0 = v0 - p;
        auto pv1 = v1 - p;
        auto pv2 = v2 - p;

        bool side1 = cross(pv0, pv1) > 0;
        bool side2 = cross(pv1, pv2) > 0;
        bool side3 = cross(pv2, pv0) > 0;

        // Point is inside if it is on the same side of all edges
        return (side1 == side2) && (side2 == side3);
    }
    static bool CalcBarycentricCoord(glm::ivec2 const& a, glm::ivec2 const& b, glm::ivec2 const& c, glm::ivec2 const& p,
                                     glm::vec3& outVal) {
        auto cross = [](glm::ivec2 const& a, glm::ivec2 const& b) { return a.x * b.y - b.x * a.y; };

        auto pa = a - p;
        auto pb = b - p;
        auto pc = c - p;

        auto s0 = cross(pa, pb);
        auto s1 = cross(pb, pc);
        auto s2 = cross(pc, pa);

        auto k0 = s0 >= 0;
        auto k1 = s1 >= 0;
        auto k2 = s2 >= 0;
        if (k0 != k1 || k1 != k2) return false;

        auto ab = b - a;
        auto ac = c - a;
        auto A = cross(ab, ac);

        outVal.x = float(s1) / A;
        outVal.y = float(s2) / A;
        outVal.z = float(s0) / A;
        return true;
    }

    void RasterizeLine(frag_in_t const& a, frag_in_t const& b) {
        auto p0 = a.fragCoord;
        auto p1 = b.fragCoord;

        auto dy = p1.y - p0.y;
        auto dx = p1.x - p0.x;

        if (std::abs(dy) > std::abs(dx)) {
            if (p0.y > p1.y)
                RasterizeLineHigh(b, a);
            else
                RasterizeLineHigh(a, b);
        } else {
            if (p0.x > p1.x)
                RasterizeLineLow(b, a);
            else
                RasterizeLineLow(a, b);
        }
    }
    void RasterizeTriangle(frag_in_t const& a, frag_in_t const& b, frag_in_t const& c) {
        auto minX = std::min(std::min(a.fragCoord.x, b.fragCoord.x), c.fragCoord.x);
        auto minY = std::min(std::min(a.fragCoord.y, b.fragCoord.y), c.fragCoord.y);
        auto maxX = std::max(std::max(a.fragCoord.x, b.fragCoord.x), c.fragCoord.x);
        auto maxY = std::max(std::max(a.fragCoord.y, b.fragCoord.y), c.fragCoord.y);

        for (auto j = minY; j <= maxY; ++j) {
            for (auto i = minX; i <= maxX; ++i) {
                glm::vec3 factor;
                if (CalcBarycentricCoord({a.fragCoord.x, a.fragCoord.y}, {b.fragCoord.x, b.fragCoord.y},
                                         {c.fragCoord.x, c.fragCoord.y}, {i, j}, factor)) {
                    frag_in_t fragIn;
                    fragIn.fragCoord.x = i;
                    fragIn.fragCoord.y = j;
                    fragIn.fragCoord.w = 1;
                    if (isPerspective) {
                        auto k = 1. / (factor.x * a.fragCoord.w + factor.y * b.fragCoord.w + factor.z * c.fragCoord.w);
                        fragIn.fragCoord.z =
                            (a.fragCoord.z * factor.x * a.fragCoord.w + b.fragCoord.z * factor.y * b.fragCoord.w +
                             c.fragCoord.z * factor.z * c.fragCoord.w) *
                            k;
                        fragIn.vertOut = (a.vertOut * factor.x * a.fragCoord.w + b.vertOut * factor.y * b.fragCoord.w +
                                          c.vertOut * factor.z * c.fragCoord.w) *
                                         k;
                    } else {
                        fragIn.fragCoord.z =
                            factor.x * a.fragCoord.z + factor.y * b.fragCoord.z + factor.z * c.fragCoord.z;
                        fragIn.vertOut = a.vertOut * factor.x + b.vertOut * factor.y + c.vertOut * factor.z;
                    }

                    Shading(fragIn);
                }
            }
        }
    }

    std::optional<per_vertex_t> ClipPoint(per_vertex_t const& a) {
        if (a.pos.x < -1 || a.pos.x > 1 || a.pos.y < -1 || a.pos.y > 1 || a.pos.z < -1 || a.pos.z > 1)
            return std::nullopt;
        return a;
    }

    // Outcode bits for homogeneous clipping planes
    enum RegionFlag: uint8_t {
        INSIDE = 0,
        LEFT = 1 << 0,
        RIGHT = 1 << 1,
        BOTTOM = 1 << 2,
        TOP = 1 << 3,
        NEAR = 1 << 4,
        FAR = 1 << 5
    };

    // Calculate outcode for a point
    uint8_t computeRegionFlags(const glm::vec4& p) {
        uint8_t code = INSIDE;

        if (p.x < -p.w)
            code |= LEFT;
        else if (p.x > p.w)
            code |= RIGHT;

        if (p.y < -p.w)
            code |= BOTTOM;
        else if (p.y > p.w)
            code |= TOP;

        if (p.z < -p.w)
            code |= NEAR;
        else if (p.z > p.w)
            code |= FAR;

        return code;
    }
    // Clip line segment against a single plane
    bool clipAgainstPlane(uint8_t plane, per_vertex_t& p0, per_vertex_t& p1) {
        float t = 0.0;
        bool intersect = false;

        // Check if line crosses the plane
        // auto needsClipping = [plane](const glm::vec4& p) { return (computeOutCode(p) & plane) != 0; };
        auto regions0 = computeRegionFlags(p0.pos);
        auto regions1 = computeRegionFlags(p1.pos);

        const bool p0Out = regions0 & plane;
        const bool p1Out = regions1 & plane;

        // Both points outside - discard
        if (p0Out && p1Out) return false;

        // Both points inside - keep
        if (!p0Out && !p1Out) return true;

        // Calculate intersection parameter t
        const glm::vec4 delta = p1.pos - p0.pos;

        if (plane == LEFT || plane == RIGHT) {
            const double sign = (plane == LEFT) ? -1.0 : 1.0;
            const double numerator = sign * p0.pos.w - p0.pos.x;
            const double denominator = delta.x - sign * delta.w;
            if (denominator != 0) t = numerator / denominator;
        } else if (plane == BOTTOM || plane == TOP) {
            const double sign = (plane == BOTTOM) ? -1.0 : 1.0;
            const double numerator = sign * p0.pos.w - p0.pos.y;
            const double denominator = delta.y - sign * delta.w;
            if (denominator != 0) t = numerator / denominator;
        } else {  // NEAR or FAR
            const double sign = (plane == NEAR) ? -1.0 : 1.0;
            const double numerator = sign * p0.pos.w - p0.pos.z;
            const double denominator = delta.z - sign * delta.w;
            if (denominator != 0) t = numerator / denominator;
        }

        if (t >= 0 && t <= 1) {
            intersect = true;
            auto intersection = p0.pos + t * delta;
            auto vertOut = p0.vertOut * (1 - t) + p1.vertOut * t;

            // Replace the outside point with intersection
            if (p0Out) {
                p0.pos = intersection;
                p0.vertOut = vertOut;
            } else {
                p1.pos = intersection;
                p1.vertOut = vertOut;
            }
        }

        return intersect;
    }

    // Nicholl–Lee–Nicholl algorithm
    bool ClipLine(per_vertex_t& a, per_vertex_t& b) {
        RegionFlag planes[]{LEFT, RIGHT, BOTTOM, TOP, NEAR, FAR};
        for (auto plane : planes) {
            if (!clipAgainstPlane(plane, a, b)) return false;
        }
        return true;
    }
    std::vector<per_vertex_t> ClipPolygonAgainstPlane(uint8_t plane, std::vector<per_vertex_t> const& vertices) {
        std::vector<per_vertex_t> ret;

        for (int i = 0; i < vertices.size(); ++i) {
            auto& p0 = vertices[i];
            auto& p1 = vertices[(i + 1) % vertices.size()];
            float t = 0.0;

            // Check if line crosses the plane
            // auto needsClipping = [plane](const glm::vec4& p) { return (computeOutCode(p) & plane) != 0; };
            auto regions0 = computeRegionFlags(p0.pos);
            auto regions1 = computeRegionFlags(p1.pos);

            const bool p0Out = regions0 & plane;
            const bool p1Out = regions1 & plane;

            // Both points outside - discard
            if (p0Out && p1Out) continue;

            // Both points inside - keep
            if (!p0Out && !p1Out) {
                // ret.emplace_back(p0);
                ret.emplace_back(p1);
                continue;
            }

            // if (!p0Out) ret.emplace_back(p0);

            // Calculate intersection parameter t
            const glm::vec4 delta = p1.pos - p0.pos;

            if (plane == LEFT || plane == RIGHT) {
                const double sign = (plane == LEFT) ? -1.0 : 1.0;
                const double numerator = sign * p0.pos.w - p0.pos.x;
                const double denominator = delta.x - sign * delta.w;
                if (denominator != 0) t = numerator / denominator;
            } else if (plane == BOTTOM || plane == TOP) {
                const double sign = (plane == BOTTOM) ? -1.0 : 1.0;
                const double numerator = sign * p0.pos.w - p0.pos.y;
                const double denominator = delta.y - sign * delta.w;
                if (denominator != 0) t = numerator / denominator;
            } else {  // NEAR or FAR
                const double sign = (plane == NEAR) ? -1.0 : 1.0;
                const double numerator = sign * p0.pos.w - p0.pos.z;
                const double denominator = delta.z - sign * delta.w;
                if (denominator != 0) t = numerator / denominator;
            }

            if (t >= 0 && t <= 1) {
                per_vertex_t p;
                p.pos = p0.pos + t * delta;
                p.vertOut = p0.vertOut * (1 - t) + p1.vertOut * t;

                ret.emplace_back(p);
                // Replace the outside point with intersection
                if (!p1Out) {
                    ret.emplace_back(p1);
                }
            }
        }
        return ret;
    }

    std::vector<per_vertex_t> ClipTrinagle(per_vertex_t const& a, per_vertex_t const& b, per_vertex_t const& c) {
        std::vector<per_vertex_t> ret;
        if (enableCullFace) {
            auto p0 = glm::vec3(Clip2NDC(a.pos));
            auto p1 = glm::vec3(Clip2NDC(b.pos));
            auto p2 = glm::vec3(Clip2NDC(c.pos));
            auto n = glm::dot(glm::cross(p1 - p0, p2 - p1), glm::vec3(0, 0, 1));
            if (n <= 0) return ret;
        }

        RegionFlag planes[]{LEFT, RIGHT, BOTTOM, TOP, NEAR, FAR};
        std::vector<per_vertex_t> polygon{a, b, c};
        for (auto plane : planes) {
            polygon = ClipPolygonAgainstPlane(plane, polygon);
            if (polygon.empty()) return {};
        }
        for (int i = 1; i < polygon.size() - 1; ++i) {
            ret.emplace_back(polygon[0]);
            ret.emplace_back(polygon[i]);
            ret.emplace_back(polygon[i + 1]);
        }
        return ret;
    }

    void DrawPoint(per_vertex_t const& a) {
        auto v = ClipPoint(a);
        if (v.has_value()) {
            // transform to screen coordinate
            frag_in_t frag;
            frag.fragCoord = this->CalcFragCoord(v->pos);
            frag.vertOut = v->vertOut;
            this->RasterizePoint(frag);
        }
    }
    void DrawLine(per_vertex_t& a, per_vertex_t& b) {
        if (ClipLine(a, b)) {
            isPerspective = a.pos.w != 1 && b.pos.w != 1;
            frag_in_t frag0;
            frag0.fragCoord = CalcFragCoord(a.pos);
            frag0.vertOut = a.vertOut;

            frag_in_t frag1;
            frag1.fragCoord = CalcFragCoord(b.pos);
            frag1.vertOut = b.vertOut;
            RasterizeLine(frag0, frag1);
        }
    }
    void DrawTriangle(per_vertex_t const& a, per_vertex_t const& b, per_vertex_t const& c) {
        auto primitives = ClipTrinagle(a, b, c);
        if (!primitives.empty()) {
            isPerspective = a.pos.w != 1 && b.pos.w != 1 && c.pos.w != 1;

            for (int i = 0; i < primitives.size(); i += 3) {
                frag_in_t frag0;
                frag0.fragCoord = CalcFragCoord(primitives[i].pos);
                frag0.vertOut = primitives[i].vertOut;

                frag_in_t frag1;
                frag1.fragCoord = CalcFragCoord(primitives[i + 1].pos);
                frag1.vertOut = primitives[i + 1].vertOut;

                frag_in_t frag2;
                frag2.fragCoord = CalcFragCoord(primitives[i + 2].pos);
                frag2.vertOut = primitives[i + 2].vertOut;
                RasterizeTriangle(frag0, frag1, frag2);
            }
        }
    }

    void VertexProcess(PrimitiveTopology topology, std::vector<per_vertex_t>& vertices) {
        // vertex post processing
        // 1.transform feedback
        // 2.flat/smooth shading
        // 3.primitive clip

        switch (topology) {
            case PRIMITIVE_TOPOLOGY_POINT_LIST:
                for (int i = 0; i < vertices.size(); ++i) DrawPoint(vertices[i]);
                break;
            case PRIMITIVE_TOPOLOGY_LINE_LIST:
                for (int i = 0; i < int(vertices.size() - 1); i += 2) DrawLine(vertices[i], vertices[i + 1]);
                break;
            case PRIMITIVE_TOPOLOGY_LINE_LOOP:
                for (int i = 0; i < int(vertices.size() - 1); ++i) DrawLine(vertices[i], vertices[i + 1]);
                break;
            case PRIMITIVE_TOPOLOGY_LINE_STRIP:
                for (int i = 0; i < int(vertices.size() - 1); ++i) DrawLine(vertices[i], vertices[i + 1]);
                DrawLine(vertices[vertices.size() - 1], vertices[0]);
                break;
            case PRIMITIVE_TOPOLOGY_TRIANGLE_LIST:
                for (int i = 0; i < int(vertices.size() - 2); i += 3)
                    DrawTriangle(vertices[i], vertices[i + 1], vertices[i + 2]);
                break;
            case PRIMITIVE_TOPOLOGY_TRIANGLE_STRIP:
                for (int i = 0; i < int(vertices.size() - 2); ++i) {
                    if (i % 2)
                        DrawTriangle(vertices[i], vertices[i + 2], vertices[i + 1]);
                    else
                        DrawTriangle(vertices[i], vertices[i + 1], vertices[i + 2]);
                }
                break;
            case PRIMITIVE_TOPOLOGY_TRIANGLE_FAN:
                for (int i = 1; i < vertices.size() - 1; ++i) DrawTriangle(vertices[0], vertices[i], vertices[i + 1]);
                break;
        }
    }

    void Draw(PrimitiveTopology topology, std::vector<vertex_attributes_t> const& verticesIn) {
        std::vector<per_vertex_t> vertices;

        // vertex process
        for (auto& v : verticesIn) vertices.emplace_back(shader.vert(v));
        VertexProcess(topology, vertices);
    }
    void DrawIndices(PrimitiveTopology topology, std::vector<vertex_attributes_t> const& verticesIn,
                     std::vector<unsigned int> const& indices) {
        std::vector<per_vertex_t> vertices;

        // vertex process
        for (auto& e : indices) {
            if (e == ~(0u)) {
                VertexProcess(topology, vertices);
                vertices.clear();
            } else
                vertices.emplace_back(shader.vert(verticesIn[e]));
        }
    }
};