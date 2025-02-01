#pragma once
#include <vector>

#include "glm/glm.hpp"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

template <typename T>
struct Image {
    std::vector<T> data;
    int width;
    int height;

    Image() = default;
    Image(int width, int height) : data(width * height), width(width), height(height) {}

    T const& operator()(int i, int j) const { return data[j * width + i]; }
    T& operator()(int i, int j) { return data[j * width + i]; }

    void fill(T const& val) { data.assign(data.size(), val); }

    void SavePNG(const char* filename) {
        stbi_flip_vertically_on_write(true);
        stbi_write_png(filename, width, height, sizeof(T), (unsigned char*)data.data(), 0);
    }

    T Sample(float x, float y) const {
        x = x - std::floor(x);
        y = y - std::floor(y);

        x *= width;
        y *= height;
        int i = x;
        int j = y;

        float cx = x - i;
        float cy = y - j;

        return (1 - cy) * ((1 - cx) * this->operator()(i, j) + cx * this->operator()(i + 1, j)) +
               cy * ((1 - cx) * this->operator()(i, j + 1) + cx * this->operator()(i + 1, j + 1));
    }
};

// struct ImageRGBA8 : Image<unsigned char>;
struct ImageRGBA8 : Image<glm::u8vec4> {
    ImageRGBA8() = default;
    ImageRGBA8(int w, int h) : Image<glm::u8vec4>(w, h) {}

    // only for u8vec4
    bool Load(const char* file) {
        int comp;
        stbi_set_flip_vertically_on_load(true);
        auto p = stbi_load(file, &width, &height, &comp, 4);
        if (p) {
            data.resize(width * height);
            memcpy(data.data(), p, data.size() * 4);
        } else {
            fprintf(stderr, "failed to load image %s\n", file);
            return false;
        }
        stbi_image_free(p);
        return true;
    }

    glm::vec4 Sample(float x, float y) const {
        x = x - std::floor(x);
        y = y - std::floor(y);

        x *= (width - 1);
        y *= (height - 1);
        int i = x;
        int j = y;

        float cx = x - i;
        float cy = y - j;

        auto getCol = [this](int i, int j) { return glm::vec4(this->operator()(i, j)) / 255.f; };
        return (1 - cy) * ((1 - cx) * getCol(i, j) + cx * getCol(std::min(i + 1, width - 1), j)) +
               cy * ((1 - cx) * getCol(i, std::min(j + 1, height - 1)) +
                     cx * getCol(i + 1, std::min(j + 1, height - 1)));
    }
};
