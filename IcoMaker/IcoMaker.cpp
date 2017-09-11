#define _WIN32_WINNT _WIN32_WINNT_WS03

#include <sdkddkver.h>

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include <windows.h>

#include <set>

#include <png.h>

#undef RGB

using namespace std;

#pragma pack(push,1)

// From the MSDN page on icon file format, there doesn't seem to be an official header for this

typedef struct
{
	BYTE        bWidth;          // Width, in pixels, of the image
	BYTE        bHeight;         // Height, in pixels, of the image
	BYTE        bColorCount;     // Number of colors in image (0 if >=8bpp)
	BYTE        bReserved;       // Reserved ( must be 0)
	WORD        wPlanes;         // Color Planes
	WORD        wBitCount;       // Bits per pixel
	DWORD       dwBytesInRes;    // How many bytes in this resource?
	DWORD       dwImageOffset;   // Where in the file is this image?
} ICONDIRENTRY, *LPICONDIRENTRY;

typedef struct
{
	WORD           idReserved;   // Reserved (must be 0)
	WORD           idType;       // Resource Type (1 for icons)
	WORD           idCount;      // How many images?
} ICONDIR, *LPICONDIR;

class RGB {
public:
	union {
		struct { BYTE b, g, r, x; };
		UINT color;
	}; 
	
	__inline RGB(int r, int g, int b)
	{
		this->r = r;
		this->g = g;
		this->b = b;
	}
	
	__inline RGB(UINT rgb)
	{
		this->color = rgb;
	}

	__inline RGB()
	{
		this->color = 0;
	}
};

class HSV 
{
public:
	union {
		struct { UINT h : 10, s : 10, v : 10, x : 2; };
		UINT color;
	};

	__inline HSV(int h, int s, int v)
	{
		this->h = h;
		this->s = s;
		this->v = v;
	}

};

#pragma pack(pop)

__inline bool operator==(const RGB& lhs, const RGB& rhs) { return lhs.color == rhs.color; }
__inline bool operator==(const HSV& lhs, const HSV& rhs) { return lhs.color == rhs.color; }
__inline bool operator<(const HSV& lhs, const HSV& rhs) { return lhs.color < rhs.color; }

class RGBFloat
{
public:
	float r;
	float g;
	float b;

	__inline RGBFloat(float r, float g, float b)
	{
		this->r = r;
		this->g = g;
		this->b = b;
	}

	__inline RGBFloat()
	{
		this->r = 0;
		this->g = 0;
		this->b = 0;
	}
};

__inline RGBFloat operator+(const RGBFloat& lhs, const RGBFloat& rhs)
{
	return RGBFloat(lhs.r + rhs.r, lhs.g + rhs.g, lhs.b + rhs.b);
}

static RGB Palette4[] = {
	0x000000,	0x800000,	0x008000,	0x808000,	0x000080,	0x800080,	0x008080,	0xC0C0C0,
	0x808080,	0xFF0000,	0x00FF00,	0xFFFF00,	0x0000FF,	0xFF00FF,	0x00FFFF,	0xFFFFFF,
};

enum class IconMode
{
	Automatic,
	Indexed16,
	Indexed256,
	Rgb24,
	Argb32,
	Png
};

class IconData
{
public:
	IconMode mode;
	int width;
	int height;
	vector<RGB> colors;
	vector<HSV> palette;
	vector<BYTE> transparency;

	vector<CHAR> raw_data;
};

enum class ResampleMethod
{
	Nearest,
	Bilinear
	// TODO: Bicubic
};

HSV rgb2hsv(RGB in)
{
	int inr = in.r;
	int ing = in.g;
	int inb = in.b;

	int minv = min(min(inr, ing), inb);
	int maxv = max(max(inr, ing), inb);

	int delta = maxv - minv;

	int outv = maxv;
	if (delta <= 0 || maxv <= 0)
		return HSV(0, 0, outv);

	int outs = delta * 256 / maxv;

	int outh;
	if (inr >= maxv)
		outh = (ing - inb) * 60 / delta;
	else if (ing >= maxv)
		outh = 120 + (inb - inr) * 60 / delta;
	else
		outh = 240 + (inr - ing) * 60 / delta;

	if (outh < 0.0)
		outh += 360.0;

	return  HSV(outh, outs, outv);
}

RGB hsv2rgb(HSV in)
{
	int inv = min(in.v, 255);
	int ins = min(in.s, 255);
	int inh = in.h;

	int outr, outg, outb;

	if (ins <= 0)
	{
		outr = inv;
		outg = inv;
		outb = inv;
		return RGB(outb, outg, outr);
	}

	if (inh >= 360) inh -= 360;

	int section = inh / 60;

	int ff = inh % 60;

	int p = inv * static_cast<int>(255 - ins) / 255;
	int q = inv * static_cast<int>(255 - ins * ff / 60) / 255;
	int t = inv * static_cast<int>(255 - ins * (60 - ff) / 60) / 255;

	switch (section) {
	case 0:
		outr = inv;
		outg = t;
		outb = p;
		break;
	case 1:
		outr = q;
		outg = inv;
		outb = p;
		break;
	case 2:
		outr = p;
		outg = inv;
		outb = t;
		break;

	case 3:
		outr = p;
		outg = q;
		outb = inv;
		break;
	case 4:
		outr = t;
		outg = p;
		outb = inv;
		break;
	case 5:
	default:
		outr = inv;
		outg = p;
		outb = q;
		break;
	}

	return RGB(outr,outg,outb);
}

static size_t indexof(const vector<RGB>& collection, RGB element)
{
	auto ptr = collection.data();
	auto size = collection.size();
	auto val = element;
	for (size_t i = 0; i < size; i++)
	{
		if (ptr[i] == val)
			return i;
	}
	return -1;
}

static size_t closestColor(const vector<HSV>& collection, RGB element)
{
	auto ptr = collection.data();
	auto size = collection.size();
	auto val = rgb2hsv(element);

	int vh = val.h;
	int vs = val.s;
	int vv = val.v;

	int diff = 0x7fffffff;
	int best = -1;
	for (size_t i = 0; i < size; i++)
	{
		if (ptr[i] == val)
			return i;

		auto cur = ptr[i];
		int ch = cur.h;
		int cs = cur.s;
		int cv = cur.v;

		int dh = vh == ch ? 0 : min(abs(vh - ch), min(abs(ch + 360 - vh), abs(vh + 360 - ch)));
		int ds = abs(vs - cs);
		int dv = abs(vv - cv);

		int df = dv * 6 + dh * 3 + ds;
		if (df < diff)
		{
			diff = df;
			best = i;
		}
	}
	return best;
}

RGB color_lerp(RGB a, RGB b, int t, int w)
{
	if (t == 0 || a == b)
		return a;

	if (t == w)
		return b;

	int ra = a.r;
	int ga = a.g;
	int ba = a.b;

	int rb = b.r;
	int gb = b.g;
	int bb = b.b;

	int ro = ra + t * (rb - ra) / w;
	int go = ga + t * (gb - ga) / w;
	int bo = ba + t * (bb - ba) / w;

	return RGB(ro, go, bo);
}

RGBFloat color_mul(RGB a, float mul)
{
	if (mul == 0)
		return RGBFloat();

	int ra = a.r;
	int ga = a.g;
	int ba = a.b;

	float ro = ra * mul;
	float go = ga * mul;
	float bo = ba * mul;

	return RGBFloat(ro, go, bo);
}

RGB color16_div(RGBFloat a, float div)
{
	if (div == 0)
		return RGB();

	float ra = a.r;
	float ga = a.g;
	float ba = a.b;

	float ro = ra / div;
	float go = ga / div;
	float bo = ba / div;

	return RGB(ro, go, bo);
}

BYTE lerp(BYTE a, BYTE b, float t)
{
	if (t < 0.0001 || a == b)
		return a;

	if (t > 0.9999)
		return b;

	int c = int(a + t * (b - a));

	return static_cast<BYTE>(c);
}

RGB color_avg(RGB a, RGB b)
{
	int ra = a.r;
	int ga = a.g;
	int ba = a.b;

	int rb = b.r;
	int gb = b.g;
	int bb = b.b;

	int ro = (rb + ra) / 2;
	int go = (gb + ga) / 2;
	int bo = (bb + ba) / 2;

	return RGB(ro, go, bo);
}

int clamp_rgb(double c)
{
	if (c < 0) return 0;
	if (c > 255) return 255;
	return static_cast<int>(c);
}

RGB swap_b_r(RGB c)
{
	return RGB(c.b,c.g,c.r);
}

static void ExtractPalette(IconData& data)
{
	data.palette.clear();

	if (data.mode == IconMode::Indexed16)
	{
		for (int i = 0; i < 16; i++)
			data.palette.push_back(rgb2hsv(Palette4[i]));
	}
	else
	{
		auto pal = set<HSV>();

		const HSV black = rgb2hsv(RGB(0));
		const HSV white = rgb2hsv(RGB(0xFFFFFF));

		pal.insert(black);
		pal.insert(white);

		for (auto argb : data.colors)
		{
			pal.insert(rgb2hsv(argb));
		}

		pal.erase(black);
		pal.erase(white);

		data.palette.push_back(black);
		data.palette.push_back(white);
		for (auto argb : pal)
			data.palette.push_back(argb);
	}
}

IconData Sharpen(IconData& source)
{
	const static int weight[9] = {
		1,  1,  1,
		1, -8,  1,
		1,  1,  1 };
	const static double alpha = 0.2;

	IconData data;
	data.width = source.width;
	data.height = source.height;
	data.mode = source.mode;
	data.colors.reserve(data.width*data.height);
	data.transparency = source.transparency;

	for (int y = 0; y < data.height; y++)
	{
		for (int x = 0; x < data.width; x++)
		{
			double pa = 0.0;
			double pr = 0.0;
			double pg = 0.0;
			double pb = 0.0;
			for (int j = -1; j < 2; j++)
			{
				for (int i = -1; i < 2; i++)
				{
					double w = weight[j * 3 + i + 4];
					int yj = max(0, min(source.height - 1, y + j));
					int xi = max(0, min(source.width - 1, x + i));
					RGB cc = source.colors[yj*data.width + xi];
					int cr = cc.r;
					int cg = cc.g;
					int cb = cc.b;
					pr += w * cr;
					pg += w * cg;
					pb += w * cb;
				}
			}

			RGB sc = source.colors[y*data.width + x];
			int sr = sc.r;
			int sg = sc.g;
			int sb = sc.b;

			int nr = clamp_rgb(sr - alpha * pr);
			int ng = clamp_rgb(sg - alpha * pg);
			int nb = clamp_rgb(sb - alpha * pb);
			data.colors.push_back(RGB(nr,ng,nb));
		}
	}

	if (data.mode == IconMode::Indexed16 || data.mode == IconMode::Indexed256)
	{
		ExtractPalette(data);
	}

	return data;
}

// Pyramid downsample (box downsample)
IconData Downsample(IconData& src, int width, int height)
{
	IconData source = src;

	while(true)
	{
		IconData data;
		data.width = source.width/2;
		data.height = source.height/2;
		data.mode = source.mode;
		data.colors.reserve(data.width*data.height);
		data.transparency.reserve(data.width*data.height);
		data.palette = source.palette;
	
		int sw = source.width - 1;
		int sh = source.height - 1;
		int tw = data.width - 1;
		int th = data.height - 1;
		for (int y = 0; y < data.height; y++)
		{
			for (int x = 0; x < data.width; x++)
			{
				int tx = x * sw;
				int ox = tx / tw;

				int ty = y * sh;
				int oy = ty / th;

				int ox2 = min(ox + 1, sw - 1);
				int oy2 = min(oy + 1, sh - 1);

				int oi0 = oy * source.width + ox;
				int oi1 = oy * source.width + ox2;
				int oi2 = oy2 * source.width + ox;
				int oi3 = oy2 * source.width + ox2;

				RGB clr0 = source.colors[oi0];
				RGB clr1 = source.colors[oi1];
				RGB clr2 = source.colors[oi2];
				RGB clr3 = source.colors[oi3];
				int a0 = source.transparency[oi0];
				int a1 = source.transparency[oi1];
				int a2 = source.transparency[oi2];
				int a3 = source.transparency[oi3];

				RGBFloat clr0a = color_mul(clr0, a0);
				RGBFloat clr1a = color_mul(clr1, a1);
				RGBFloat clr2a = color_mul(clr2, a2);
				RGBFloat clr3a = color_mul(clr3, a3);

				float cct = a0 + a1 + a2 + a3;
				RGB vc = color16_div(clr0a + clr1a + clr2a + clr3a, cct);
				data.colors.push_back(vc);

				int va = cct/4;
				data.transparency.push_back(va);
			}
		}

		// Not needed the last loop
		if (data.width <= width)
			return data;

		source = data;
	}
}

int is_power_of_two(UINT n)
{
	UINT a = ((n & 0xAAAAAAAA) >> 1) + (n & 0x55555555);
	UINT b = ((a & 0xCCCCCCCC) >> 2) + (a & 0x33333333);
	UINT c = ((b & 0xF0F0F0F0) >> 4) + (b & 0x0F0F0F0F);
	UINT d = ((c & 0xFF00FF00) >> 8) + (c & 0x00FF00FF);
	UINT e = ((d & 0xFFFF0000) >> 16) + (d & 0x0000FFFF);
	return e == 1;
}

int next_power_of_two(UINT n)
{
	n = n | (n >> 1);
	n = n | (n >> 2);
	n = n | (n >> 4);
	n = n | (n >> 8);
	n = n | (n >> 16);
	return n+1;
}

IconData Resample(IconData& src, int width, int height, ResampleMethod method = ResampleMethod::Bilinear)
{
	IconData source;

	// TODO: Vertical-only and Horizontal-only Downsampling, and downsample before 
	if (src.width == src.height && is_power_of_two(src.width) && width < src.width && height < src.height)
	{
		if (width == height && is_power_of_two(width))
		{
			return Downsample(src, width, height);
		}

		int temp_w = next_power_of_two(width);
		int temp_h = next_power_of_two(height);
		if (temp_w < src.width || temp_h < src.height)
		{
			source = Downsample(src, temp_w, temp_h);
		}
		else
		{
			source = src;
		}
	}
	else
	{
		source = src;
	}


	// Pre-sharpening step, helps avoid some of the incuded blur
	//if (source.mode == IconMode::Indexed16)
	//	ExtractPalette(source);

	IconData data;
	data.width = width;
	data.height = height;
	data.mode = source.mode;
	data.colors.reserve(width*height);
	data.transparency.reserve(width*height);
	data.palette = source.palette;

	int offx = width / 2;
	int offy = height / 2;
	if (method == ResampleMethod::Nearest)
	{
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				int tx = x * source.width;
				int ty = y * source.height;
				int ox = (tx + offx) / width;
				int oy = (ty + offy) / height;
				int oi = oy * source.width + ox;
				data.colors.push_back(source.colors[oi]);
				data.transparency.push_back(source.transparency[oi]);
			}
		}
	}
	else
	{
		int sw = source.width - 1;
		int sh = source.height - 1;
		int tw = width - 1;
		int th = height - 1;
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				int tx = x * sw;
				int ox = tx / tw;
				float fx = tx/float(tw) - ox;

				int ty = y * sh;
				int oy = ty / th;
				float fy = ty / float(th) - oy;

				int ox2 = min(ox + 1, sw - 1);
				int oy2 = min(oy + 1, sh - 1);

				int oi0 = oy * source.width + ox;
				int oi1 = oy * source.width + ox2;
				int oi2 = oy2 * source.width + ox;
				int oi3 = oy2 * source.width + ox2;

				RGB clr0 = source.colors[oi0];
				RGB clr1 = source.colors[oi1];
				RGB clr2 = source.colors[oi2];
				RGB clr3 = source.colors[oi3];
				int a0 = source.transparency[oi0];
				int a1 = source.transparency[oi1];
				int a2 = source.transparency[oi2];
				int a3 = source.transparency[oi3];

				float cc0 = a0 * (1 - fx) * (1 - fy);
				float cc1 = a1 * fx * (1 - fy);
				float cc2 = a2 * (1 - fx) * fy;
				float cc3 = a3 * fx * fy;
				float cct = cc0+cc1+cc2+cc3;

				RGBFloat clr0a = color_mul(clr0, cc0);
				RGBFloat clr1a = color_mul(clr1, cc1);
				RGBFloat clr2a = color_mul(clr2, cc2);
				RGBFloat clr3a = color_mul(clr3, cc3);

				RGB vc = color16_div(clr0a + clr1a + clr2a + clr3a, cct);
				data.colors.push_back(vc);

				int ha1 = a0 * (1 - fx) + a1 * fx;
				int ha2 = a2 * (1 - fx) + a3 * fx;
				int va = ha1 * (1 - fy) + ha2 * fy;
				data.transparency.push_back(va);
			}
		}
	}

	return data;
}

static void FillPalette(IconData& data)
{
	if (data.mode == IconMode::Indexed16)
	{
		for (int i = data.palette.size(); i < 16; i++)
			data.palette.push_back(rgb2hsv(Palette4[i]));
	}
	else
	{
		for (int i = data.palette.size(); i < 256; i++)
			data.palette.push_back(rgb2hsv(RGB((i * 0xFFFFFF) / 255)));
	}
}

static bool repalette_warning_shown = false;
IconData Repalette(IconData& source)
{
	IconData data;
	data.width = source.width;
	data.height = source.height;
	data.mode = source.mode;
	data.colors.reserve(data.width*data.height);
	data.transparency = source.transparency;

	if (!repalette_warning_shown)
	{
		repalette_warning_shown = true;
		cout << L"Warning: The input image had too many colors, and its palette had to be reduced. This will result in lower quality than pre-processing in a dedicated program!" << endl;
	}

#if 0
	// BROKEN: The palette is now in HSV form
	// EXTREMELY INEFFICIENT CODE LIES HERE!
	data.palette = source.palette;
	while (data.palette.size() > 256)
	{
		int ai = -1;
		RGB a = 0;
		RGB b = 0;
		int d = 1 << 24;
		for (int y = 0; y < data.palette.size(); y++)
		{
			RGB cy = data.palette[y];
			for (int x = y + 1; x < data.palette.size(); x++)
			{
				RGB cx = data.palette[x];

				int df = colordiff(cy, cx);
				if (df < d)
				{
					d = df;
					a = cy;
					b = cx;
					ai = y;
				}
			}
		}
		data.palette[ai] = color_avg(a, b);
		data.palette.erase(data.palette.end() - 1);
	}
#else
	int forced = 2;
	int nsource = source.palette.size() - 1 - forced;
	int ntotal = 255 - forced;

	for (int i = 0; i < forced; i++)
	{
		data.palette.push_back(source.palette[i]);
	}
	for (int i = forced; i < 256; i++)
	{
		data.palette.push_back(source.palette[forced + (i - forced) * nsource / ntotal]);
	}
#endif

	for (auto color : source.colors)
	{
		auto idx = closestColor(data.palette, color);
		data.colors.push_back(hsv2rgb(data.palette[idx]));
	}

	return data;
}

IconData ReadPngFile(wstring filename, int width = -1, int height = -1, IconMode mode = IconMode::Automatic)
{
	FILE *fp;

	if (_wfopen_s(&fp, filename.c_str(), L"rb")) abort();

	png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
	if (!png) abort();

	png_infop info = png_create_info_struct(png);
	if (!info) abort();

	if (setjmp(png_jmpbuf(png)))
		abort();

	png_init_io(png, fp);

	png_read_info(png, info);

	IconData data;
	data.mode = mode;
	data.width = png_get_image_width(png, info);
	data.height = png_get_image_height(png, info);
	auto color_type = png_get_color_type(png, info);
	auto bit_depth = png_get_bit_depth(png, info);

	// Read any color_type into 8bit depth, RGBA format.
	// See http://www.libpng.org/pub/png/libpng-manual.txt

	if (bit_depth == 16)
		png_set_strip_16(png);

	if (color_type == PNG_COLOR_TYPE_PALETTE)
	{
		if (data.mode == IconMode::Automatic)
			data.mode = IconMode::Indexed256;
		png_set_palette_to_rgb(png);
	}
	// PNG_COLOR_TYPE_GRAY_ALPHA is always 8 or 16bit depth.
	else if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
	{
		if (data.mode == IconMode::Automatic)
			data.mode = IconMode::Indexed256;
		png_set_expand_gray_1_2_4_to_8(png);
	}
	else
	{
		if (data.mode == IconMode::Automatic)
			data.mode = data.width > 64 ? IconMode::Png : IconMode::Argb32;
	}

	if (png_get_valid(png, info, PNG_INFO_tRNS))
		png_set_tRNS_to_alpha(png);

	// These color_type don't have an alpha channel then fill it with 0xff.
	if (color_type == PNG_COLOR_TYPE_RGB ||
		color_type == PNG_COLOR_TYPE_GRAY ||
		color_type == PNG_COLOR_TYPE_PALETTE)
		png_set_filler(png, 0xFF, PNG_FILLER_AFTER);

	if (color_type == PNG_COLOR_TYPE_GRAY ||
		color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
		png_set_gray_to_rgb(png);

	png_read_update_info(png, info);

	auto row_pointers = new png_bytep[data.height];
	for (int y = 0; y < data.height; y++) {
		row_pointers[y] = new png_byte[png_get_rowbytes(png, info)];
	}

	png_read_image(png, row_pointers);

	for (int y = 0; y < data.height; y++)
	{
		auto cr = reinterpret_cast<UINT*>(row_pointers[data.height - y - 1]);
		for (int x = 0; x < data.width; x++)
		{
			auto cc = cr[x];
			data.transparency.push_back(cc >> 24);
			data.colors.push_back(swap_b_r(RGB(cc)));
		}
	}

	for (int y = 0; y < data.height; y++) {
		delete[] row_pointers[y];
		row_pointers[y] = nullptr;
	}
	delete[] row_pointers;

	fclose(fp);

	if (data.mode == IconMode::Indexed16 || data.mode == IconMode::Indexed256)
	{
		ExtractPalette(data);
	}

	if ((width > 0 && width != data.width) || (height > 0 && height != data.height))
	{
		data = Resample(data, width, height);
	}

	if (data.mode == IconMode::Indexed256 && data.palette.size() > 256)
	{
		data = Repalette(data);
	}

	if (data.mode == IconMode::Indexed256)
	{
		FillPalette(data);
	}

	return data;
}

void WriteToVector(png_structp png_ptr, png_bytep pdata, png_size_t length)
{
	auto& data = *static_cast<IconData*>(png_get_io_ptr(png_ptr));
	int original_size = data.raw_data.size();
	data.raw_data.resize(original_size + length);
	memcpy(data.raw_data.data() + original_size, pdata, length);
}

void WritePngFile(IconData& data) {

	png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (!png) abort();

	png_infop info = png_create_info_struct(png);
	if (!info) abort();

	if (setjmp(png_jmpbuf(png)))
		abort();

	//png_init_io(png, fp);

	png_set_write_fn(png, &data, WriteToVector, nullptr);

	// Output is 8bit depth, RGBA format.
	png_set_IHDR(
		png,
		info,
		data.width, data.height,
		8,
		PNG_COLOR_TYPE_RGBA,
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT
	);
	png_write_info(png, info);

	auto buffRGB = data.colors.data();
	auto buffA = data.transparency.data();

	auto row_pointers = new png_bytep[data.height];
	for (int y = 0; y < data.height; y++) {
		row_pointers[y] = new png_byte[data.width * 4];
		auto dst = reinterpret_cast<UINT*>(row_pointers[y]);
		auto ptrRGB = (buffRGB + data.width * (data.height - y - 1));
		auto ptrA = (buffA + data.width * (data.height - y - 1));
		for (int i = 0; i < data.width; i++)
			dst[i] = (ptrA[i] << 24) | (ptrRGB[i].b << 16) | (ptrRGB[i].g << 8) | ptrRGB[i].r;
	}

	png_write_image(png, row_pointers);
	png_write_end(png, nullptr);

	for (int y = 0; y < data.height; y++) {
		delete[] row_pointers[y];
		row_pointers[y] = nullptr;
	}
	delete[] row_pointers;
}

int Bpp(IconMode mode)
{
	switch (mode)
	{
	case IconMode::Indexed16: return 4;
	case IconMode::Indexed256: return 8;
	case IconMode::Rgb24: return 24;
	case IconMode::Argb32: return 32;
	case IconMode::Png: return 32;
	}
	return 32;
}

bool WriteBmp(IconData& data)
{
	const bool useAndMask = true; // required for (bitsPerPixel < 32), but for compatibility, it's best to always have one.

	vector<CHAR>& outbuffer(data.raw_data);
	const vector<RGB>& pixels(data.colors);
	const int width = data.width;
	const int height = data.height;
	const WORD bitsPerPixel = Bpp(data.mode);
	const vector<HSV>& palette(data.palette);
	const vector<BYTE>& transparency(data.transparency);

	DWORD palLength = palette.size() == (1 << bitsPerPixel) ? 0 : palette.size();

	DWORD palSize = palette.size() * 4;
	DWORD pixelsSizeXor = ((width * bitsPerPixel + 31) / 32) * 4 * height;
	DWORD pixelsSizeAnd = useAndMask ? (((width + 31) / 32) * 4 * height) : 0;
	DWORD pixelsSize = pixelsSizeXor + pixelsSizeAnd;

	DWORD palOffset = sizeof(BITMAPINFOHEADER);
	DWORD pixelsOffset = palOffset + palSize;
	DWORD pixelsOffsetAnd = pixelsOffset + pixelsSizeXor;
	DWORD totalSize = pixelsOffset + pixelsSize;

	char* buffer = new char[totalSize];

	*reinterpret_cast<BITMAPINFOHEADER*>(buffer) = {
		40,
		width,
		useAndMask ? height * 2 : height,
		1,
		bitsPerPixel,
		0,
		pixelsSize,
		0,
		0,
		palLength,
		palLength
	};

	if (palSize > 0)
	{
		auto bb = reinterpret_cast<RGB*>(buffer + palOffset);
		for (int i = 0; i < palette.size(); i++)
			bb[i] = hsv2rgb(palette[i]);
	}

	if (bitsPerPixel == 1)
	{
		bool flush = (width & 3) != 0;
		for (int y = 0, i = 0, o = 0; y < height; y++)
		{
			int c = 0;
			int b = 0;
			for (int x = 0; x < width; x++, i++)
			{
				auto color = pixels[i];
				auto transparent = (~transparency[i] >> 7) & 1;
				auto index = transparent ? 0 : closestColor(palette, color);
				switch (x & 7)
				{
				case 0: c = index; break;
				case 1: c = (c << 1) | index; break;
				case 2: c = (c << 1) | index; break;
				case 3: c = (c << 1) | index; break;
				case 4: c = (c << 1) | index; break;
				case 5: c = (c << 1) | index; break;
				case 6: c = (c << 1) | index; break;
				default:
					c = (c << 1) | index;
					buffer[pixelsOffset + o++] = c;
					b++;
					break;
				}
			}
			if (flush)
			{
				buffer[pixelsOffset + o++] = c;
				b++;
			}
			while ((b++ % 4) != 0)
				buffer[pixelsOffsetAnd + o++] = 0;
		}
	}
	else if (bitsPerPixel == 4)
	{
		bool flush = (width & 1) != 0;
		for (int y = 0, i = 0, o = 0; y < height; y++)
		{
			int c = 0;
			int b = 0;
			for (int x = 0; x < width; x++, i++)
			{
				auto color = pixels[i];
				auto transparent = (~transparency[i] >> 7) & 1;
				auto index = transparent ? 0 : closestColor(palette, color);
				if ((x & 1) == 0)
				{
					c = index << 4;
				}
				else
				{
					c = c | index;
					buffer[pixelsOffset + o++] = c;
					b++;
				}
			}
			if (flush)
			{
				buffer[pixelsOffset + o++] = c;
				b++;
			}
			while ((b++ % 4) != 0)
				buffer[pixelsOffsetAnd + o++] = 0;
		}
	}
	else if (bitsPerPixel == 8)
	{
		for (int y = 0, i = 0, o = 0; y < height; y++)
		{
			int x;
			for (x = 0; x < width; x++, i++)
			{
				auto color = pixels[i];
				auto transparent = (~transparency[i] >> 7) & 1;
				auto index = transparent ? 0 : closestColor(palette, color);
				buffer[pixelsOffset + o++] = index;
			}
			while ((x++ % 4) != 0)
				buffer[pixelsOffsetAnd + o++] = 0;
		}
	}
	else if (bitsPerPixel == 24)
	{
		for (int y = 0, i = 0, o = 0; y < height; y++)
		{
			int x;
			for (x = 0; x < width; x++, i++)
			{
				auto color = pixels[i];
				auto transparent = (~transparency[i] >> 7) & 1;
				auto index = transparent ? 0 : color.color;
				buffer[pixelsOffset + o++] = index;
				buffer[pixelsOffset + o++] = index >> 8;
				buffer[pixelsOffset + o++] = index >> 16;
			}
			while ((x++ % 4) != 0)
				buffer[pixelsOffsetAnd + o++] = 0;
		}
	}
	else if (bitsPerPixel == 32)
	{
		auto bptr = reinterpret_cast<UINT*>(buffer + pixelsOffset);
		for (int i = 0; i < pixels.size(); i++)
		{
			bptr[i] = pixels[i].color | (transparency[i] << 24);
		}
	}
	else
	{
		wcerr << "Invalid pixel format, no data written!" << endl;
		return false;
	}

	if (useAndMask) // always true
	{
		bool flush = (width & 3) != 0;
		for (int y = 0, i = 0, o = 0; y < height; y++)
		{
			int c = 0;
			int b = 0;
			for (int x = 0; x < width; x++, i++)
			{
				auto transparent = (~transparency[i] >> 7) & 1;
				switch (x & 7)
				{
				case 0: c = transparent; break;
				case 1: c = (c << 1) | transparent; break;
				case 2: c = (c << 1) | transparent; break;
				case 3: c = (c << 1) | transparent; break;
				case 4: c = (c << 1) | transparent; break;
				case 5: c = (c << 1) | transparent; break;
				case 6: c = (c << 1) | transparent; break;
				default:
					c = (c << 1) | transparent;
					buffer[pixelsOffsetAnd + o++] = c;
					b++;
					break;
				}
			}
			if (flush)
			{
				buffer[pixelsOffsetAnd + o++] = c;
				b++;
			}
			while ((b++ % 4) != 0)
				buffer[pixelsOffsetAnd + o++] = 0;
		}
	}

	outbuffer.resize(totalSize);
	RtlCopyMemory(outbuffer.data(), buffer, totalSize);

	return true;
}

void EncodeImageData(IconData& data)
{
	if (data.mode == IconMode::Png)
	{
		//EncodePNG(data);
		WritePngFile(data);
	}
	else
	{
		WriteBmp(data);
	}
}

void WriteIcon(ofstream& outfile, vector<IconData>& images, bool isCursor = false, int hotspotX = 0, int hotspotY = 0)
{
	ICONDIR header = {
		static_cast<WORD>(0),
		static_cast<WORD>(isCursor ? 2 : 1),
		static_cast<WORD>(images.size())
	};

	outfile.write(reinterpret_cast<PCHAR>(&header), sizeof(header));

	auto startOfData = sizeof(header) + images.size() * sizeof(ICONDIRENTRY);
	auto currentImageStart = startOfData;

	for (auto& d : images)
	{
		EncodeImageData(d);

		ICONDIRENTRY entry = {
			static_cast<BYTE>(d.width & 0xFF),
			static_cast<BYTE>(d.height & 0xFF),
			static_cast<BYTE>(Bpp(d.mode) >= 8 ? 0 : d.palette.size()),
			static_cast<BYTE>(0),
			static_cast<WORD>(isCursor ? hotspotX : (d.mode == IconMode::Indexed16 ? 0 : 1)), /*planes*/
			static_cast<WORD>(isCursor ? hotspotY : (d.mode == IconMode::Indexed16 ? 0 : Bpp(d.mode))),
			static_cast<DWORD>(d.raw_data.size()),
			static_cast<DWORD>(currentImageStart)
		};

		outfile.write(reinterpret_cast<PCHAR>(&entry), sizeof(entry));

		currentImageStart += entry.dwBytesInRes;
	}

	for (auto& d : images)
	{
		outfile.write(d.raw_data.data(), d.raw_data.size());
	}
}

void WriteOutput(vector<IconData>& images, wstring output_filename)
{
	auto outfile = ofstream(output_filename.c_str(), ios::out | ios::binary);

	WriteIcon(outfile, images);

	outfile.close();
}

int wmain(int argc, wchar_t** argv)
{
	vector<IconData> datas;

	wstring input_filename;
	int input_width;
	int input_height;

	wstring output_filename;

	int state = 0;
	for (int i = 1; i < argc; )
	{
		auto arg = argv[i];

		switch (state)
		{
		case 0:
			if (arg[0] == L'-')
			{
				if (wcscmp(arg, L"-a") == 0 || wcscmp(arg, L"--add") == 0)
				{
					input_filename = L"";
					input_width = -1;
					input_height = -1;
					state = 2;
				}
				else if (wcscmp(arg, L"-o") == 0 || wcscmp(arg, L"--output") == 0)
				{
					output_filename = L"";
					state = 1;
				}
				else
				{
					wcout << L"Unknown flag: " << arg << endl;
					goto show_help;
				}
			}
			else
			{
				wcout << L"Unexpected value" << endl;
				goto show_help;
			}
			break;
		case 1:
			if (arg[0] == L'-')
			{
				wcout << L"Missing filename argument for --output" << endl;
				goto show_help;
			}

			output_filename = arg;
			state = 0;
			break;
		case 2:
			if (arg[0] == L'-')
			{
				wcout << L"Missing filename argument for --add" << endl;
				goto show_help;
			}

			input_filename = arg;
			state = 3;
			break;
		case 3:
			if (arg[0] == L'-')
			{
				datas.push_back(ReadPngFile(input_filename));
				state = 0;
				continue; // without consuming
			}

			input_width = _wtoi(arg);
			state = 4;
			break;
		case 4:
			if (arg[0] == L'-')
			{
				datas.push_back(ReadPngFile(input_filename, input_width, input_width));
				state = 0;
				continue; // without consuming
			}

			input_height = _wtoi(arg);
			state = 5;
			break;
		case 5:
			if (arg[0] == L'-')
			{
				datas.push_back(ReadPngFile(input_filename, input_width, input_height));
				state = 0;
				continue; // without consuming
			}

			IconMode input_format;
			if (wcscmp(arg, L"i4") == 0 || wcscmp(arg, L"indexed4") == 0)
				input_format = IconMode::Indexed16;
			else if (wcscmp(arg, L"i8") == 0 || wcscmp(arg, L"indexed8") == 0)
				input_format = IconMode::Indexed256;
			else if (wcscmp(arg, L"rgb") == 0)
				input_format = IconMode::Rgb24;
			else if (wcscmp(arg, L"argb") == 0)
				input_format = IconMode::Argb32;
			else if (wcscmp(arg, L"png") == 0)
				input_format = IconMode::Png;
			else
			{
				wcout << L"Unknown icon pixel format: " << arg << endl;
				goto show_help;
			}

			datas.push_back(ReadPngFile(input_filename, input_width, input_height, input_format));
			state = 0;
			break;
		}

		i++;
	}

	switch (state)
	{
	case 0:
		break;
	case 1:
		wcout << L"Missing filename argument for --output" << endl;
		goto show_help;
	case 2:
		wcout << L"Missing filename argument for --add" << endl;
		goto show_help;
	case 3:
		datas.push_back(ReadPngFile(input_filename));
		break;
	case 4:
		datas.push_back(ReadPngFile(input_filename, input_width, input_width));
		break;
	case 5:
		datas.push_back(ReadPngFile(input_filename, input_width, input_height));
		break;
	default:
		wcout << L"Unexpected state" << endl;
		goto show_help;
	}

	if (output_filename.size() == 0)
	{
		wcout << L"Output filename not provided" << endl;
		goto show_help;
	}

	if (datas.size() == 0)
	{
		wcout << L"No images provided" << endl;
		goto show_help;
	}

	WriteOutput(datas, output_filename);

	return 0;

show_help:
	wcout << L"TODO: Help lines" << endl;

	return 1;
}
