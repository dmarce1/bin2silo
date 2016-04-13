/*
 * bin2silo.cpp
 *
 *  Created on: Apr 13, 2016
 *      Author: dmarce1
 */

#include <silo.h>
#include <array>
#include <vector>
#include <stdlib.h>
#include <cstdint>
#include <climits>
#include <math.h>

/*******************************************************************/
// copied from http://stackoverflow.com/questions/105252/how-do-i-convert-between-big-endian-and-little-endian-values-in-c
template<typename T>
T swap_endian(T u) {
	static_assert (CHAR_BIT == 8, "CHAR_BIT != 8");

	union {
		T u;
		unsigned char u8[sizeof(T)];
	} source, dest;

	source.u = u;

	for (size_t k = 0; k < sizeof(T); k++)
		dest.u8[k] = source.u8[sizeof(T) - k - 1];

	return dest.u;
}
/***************************************************************/

constexpr int nf = 14;

const char* field_names[] = { "den", "frac1", "frac2", "pot", "pres", "spec1",
		"spec2", "spec3", "spec4", "spec5", "tau", "velphi", "velr", "velz" };

int main(int argc, char* argv[]) {
	int nphi, nz, nr, fnum;
	if (argc != 5) {
		printf("bin2silo <input filenumber> nphi nz nr\n");
		abort();
	}
	fnum = atoi(argv[1]);
	nphi = atoi(argv[2]);
	nz = atoi(argv[3]);
	nr = atoi(argv[4]);
	const int total_sz = nphi * nz * (nr-1);
	std::array<std::vector<double>, nf> data;
	std::uint32_t header;
	for (int f = 0; f != nf; ++f) {
		data[f].resize(total_sz);
		char* filename;
		asprintf(&filename, "%s_%i", field_names[f], fnum);
		printf("Reading %s\n", filename);
		FILE* fp = fopen(filename, "rb");
		if (fp == NULL) {
			printf("unable to open %s\n", filename);
			abort();
		}
		fread(&header, 4, 1, fp);
		header = swap_endian(header);
		for (int i = 0; i != nphi; ++i) {
			for (int j = 0; j != nz; ++j) {
				for (int k = 0; k != nr; ++k) {
					double v;
					fread(&v, sizeof(v), 1, fp);
					v = swap_endian(v);
					const int iii = i * nz * (nr-1) + j * (nr-1) + (k-1);
					if( k != 0 ) {
						data[f][iii] = v;
					}
				}
			}
		}
		free(filename);
		fclose(fp);
	}

	struct point {
		double x, y, z;
	};
	std::vector<point> pts;
	pts.reserve(nr * (nphi + 1) * (nz + 1));
	double dx = 2.0 * M_PI / double(nphi);
	std::vector<int> zone_nodes;
	zone_nodes.reserve(8 * total_sz);

	for (int i = 0; i != nphi + 1; ++i) {
		for (int j = 0; j != nz + 1; ++j) {
			for (int k = 0; k != nr; ++k) {
				double phi = (i) * dx;
				double z = (j-nz/2) * dx;
				double r = k * dx;
				point pt;
				pt.x = r * cos(phi);
				pt.y = r * sin(phi);
				//	if( r > 1.0 )
				//	printf( "%e\n", pt.y);
				pt.z = z;
				pts.push_back(pt);
			}
		}
	}
	printf("%i %i %i\n", nphi, nz, nr);
	for (int i = 0; i != nphi; ++i) {
		for (int j = 0; j != nz; ++j) {
			for (int k = 0; k != nr-1; ++k) {
				const int i0 = (i + 0) * (nr * (nz + 1))
						+ (j + 0) * nr + (k + 0);
				const int i1 = (i + 1) * (nr * (nz + 1))
						+ (j + 0) * nr + (k + 0);
				const int i2 = (i + 1) * (nr * (nz + 1))
						+ (j + 1) * nr + (k + 0);
				const int i3 = (i + 0) * (nr * (nz + 1))
						+ (j + 1) * nr + (k + 0);
				const int i4 = (i + 0) * (nr * (nz + 1))
						+ (j + 0) * nr + (k + 1);
				const int i5 = (i + 1) * (nr * (nz + 1))
						+ (j + 0) * nr + (k + 1);
				const int i6 = (i + 1) * (nr * (nz + 1))
						+ (j + 1) * nr + (k + 1);
				const int i7 = (i + 0) * (nr * (nz + 1))
						+ (j + 1) * nr + (k + 1);
				zone_nodes.push_back(i0);
				zone_nodes.push_back(i1);
				zone_nodes.push_back(i2);
				zone_nodes.push_back(i3);
				zone_nodes.push_back(i4);
				zone_nodes.push_back(i5);
				zone_nodes.push_back(i6);
				zone_nodes.push_back(i7);

			}
		}
	}
	std::vector<double> xs(pts.size());
	std::vector<double> ys(pts.size());
	std::vector<double> zs(pts.size());
	for (int i = 0; i != pts.size(); ++i) {
		xs[i] = pts[i].x;
		ys[i] = pts[i].y;
		zs[i] = pts[i].z;
	}
	DBfile *db = DBCreateReal("output.silo", DB_CLOBBER, DB_LOCAL, "Euler Mesh",
	DB_PDB);
	int nshapes = 1;
	int shapesize[1] = { 8 };
	int shapetype[1] = { DB_ZONETYPE_HEX };
	int shapecnt[1] = { total_sz };
	const char* coord_names[3] = { "x", "y", "z" };
	std::array<double*, 3> node_coords = { xs.data(), ys.data(), zs.data() };
	DBPutZonelist2(db, "zones", total_sz, 3, zone_nodes.data(),
			zone_nodes.size(), 0, 0, 0, shapetype, shapesize, shapecnt, nshapes,
			NULL);
	DBPutUcdmesh(db, "mesh", 3, coord_names, node_coords.data(), pts.size(),
			total_sz, "zones", nullptr, DB_DOUBLE, NULL);
	for (int f = 0; f != nf; ++f) {
		DBPutUcdvar1(db, field_names[f], "mesh", data[f].data(), total_sz, NULL,
				0, DB_DOUBLE, DB_ZONECENT, NULL);
	}

	DBClose(db);
}
