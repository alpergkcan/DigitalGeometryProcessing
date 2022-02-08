#pragma once
#include "pch.h"

typedef class QH {
public:
	QH(){};

	static Quat to_quat(const Mat3& mat) {
		float sq = (1 + mat.coeff(0, 0) + mat.coeff(1, 1) + mat.coeff(2, 2)) / 4;
		float w, x, y, z;
		 if (sq > FLT_EPSILON) {
			 w = sqrt(sq);
			 x = (mat.coeff(1, 2) - mat.coeff(2, 1)) / (4 * w);
			 y = (mat.coeff(2, 0) - mat.coeff(0, 2)) / (4 * w);
			 z = (mat.coeff(0, 1) - mat.coeff(1, 0)) / (4 * w);
		} else {
			w = 0;
			sq = -1 * (mat.coeff(1, 1) + mat.coeff(2, 2)) / 2;
			if (sq > FLT_EPSILON) {
				x = sqrt(sq);
				y = mat.coeff(0, 1) / (2 * x);
				z = mat.coeff(0, 2) / (2 * x);
			}else{
				x = 0;
				sq = 0.5f * (1 - mat.coeff(2, 2));
				if (sq > FLT_EPSILON) {
					y = sqrt(sq);
					z = mat.coeff(1, 2) / (2 * y);
				} else {
					y = 0;
					z = 1;
				}
			}
		}

		Quat res;
		res.w() = w;
		res.x() = x;
		res.y() = y;
		res.z() = z;
		return res.normalized();
	}

	static Mat3 to_mat(const Quat& qua) {
		//if (!(1 - FLT_EPSILON <= qua.norm() <= 1 + FLT_EPSILON)) {
		//	qua.normalize();
		//}
		
		float xx, xy, xz, xw, yy, yz, yw, zz, zw;

		xy = 2 * qua.x() * qua.y();
		xz = 2 * qua.x() * qua.z();
		xw = 2 * qua.x() * qua.w();

		yz = 2 * qua.y() * qua.z();
		yw = 2 * qua.y() * qua.w();

		zw = 2 * qua.z() * qua.w();

		xx = 2 * qua.x() * qua.x();
		yy = 2 * qua.y() * qua.y();
		zz = 2 * qua.z() * qua.z();

		Mat3 A;
		A << 1 - yy - zz,     xy + zw,     xz - yw,
			     xy - zw, 1 - xx - zz,     yz + xw,
			     xz + yw,     yz - xw, 1 - xx - yy;
		return A;
	}

	static float ReciprocalSqrt(float x) {
		long i;
		float y, r;
		y = x * 0.5f;
		i = *(long*)(&x);
		i = 0x5f3759df - (i >> 1);
		r = *(float*)(&i);
		r = r * (1.5f - r * r * y);
		return r;
	}
	static Quat to_quat_recip(Mat3 m) {

		Quat q;
		if (m.coeff(0 , 0) + m.coeff(1 , 1) + m.coeff(2 , 2) > 0.0f) {
			float t = +m.coeff(0 , 0) + m.coeff(1 , 1) + m.coeff(2 , 2) + 1.0f;
			float s = ReciprocalSqrt(t) * 0.5f;
			q.w() = s * t;
			q.z() = (m.coeff(0 , 1) - m.coeff(1 , 0)) * s;
			q.y() = (m.coeff(2 , 0) - m.coeff(0 , 2)) * s;
			q.x() = (m.coeff(1 , 2) - m.coeff(2 , 1)) * s;
		}
		else if (m.coeff(0 , 0) > m.coeff(1 , 1) && m.coeff(0 , 0) > m.coeff(2 , 2)) {
			float t = +m.coeff(0 , 0) - m.coeff(1 , 1) - m.coeff(2 , 2) + 1.0f;
			float s = ReciprocalSqrt(t) * 0.5f;
			q.x() = s * t;
			q.y() = (m.coeff(0 , 1) + m.coeff(1 , 0)) * s;
			q.z() = (m.coeff(2 , 0) + m.coeff(0 , 2)) * s;
			q.w() = (m.coeff(1 , 2) - m.coeff(2 , 1)) * s;
		}
		else if (m.coeff(1 , 1) > m.coeff(2 , 2)) {
			float t = -m.coeff(0 , 0) + m.coeff(1 , 1) - m.coeff(2 , 2) + 1.0f;
			float s = ReciprocalSqrt(t) * 0.5f;
			q.y() = s * t;
			q.x() = (m.coeff(0 , 1) + m.coeff(1 , 0)) * s;
			q.w() = (m.coeff(2 , 0) - m.coeff(0 , 2)) * s;
			q.z() = (m.coeff(1 , 2) + m.coeff(2 , 1)) * s;
		}
		else {
			float t = -m.coeff(0 , 0) - m.coeff(1 , 1) + m.coeff(2 , 2) + 1.0f;
			float s = ReciprocalSqrt(t) * 0.5f;
			q.z() = s * t;
			q.w() = (m.coeff(0 , 1) - m.coeff(1 , 0)) * s;
			q.x() = (m.coeff(2 , 0) + m.coeff(0 , 2)) * s;
			q.y() = (m.coeff(1 , 2) + m.coeff(2 , 1)) * s;
		}
		return q;
	}

	static Quat slerp(Quat& qa, Quat& qb, float t) {
		Quat res;
		float cosHalfTheta = qa.w() * qb.w() + qa.x() * qb.x() + qa.y() * qb.y() + qa.z() * qb.z();
		if (abs(cosHalfTheta) >= 1.0) {
			res = qa;
			return res;
		}
		float halfTheta = acos(cosHalfTheta);
		float sinHalfTheta = sqrt(1.0 - cosHalfTheta * cosHalfTheta);

		if (fabs(sinHalfTheta) < 0.001) { // fabs is floating point absolute
			res.w() = (qa.w() * 0.5 + qb.w() * 0.5);
			res.x() = (qa.x() * 0.5 + qb.x() * 0.5);
			res.y() = (qa.y() * 0.5 + qb.y() * 0.5);
			res.z() = (qa.z() * 0.5 + qb.z() * 0.5);
			return res;
		}
		float  ratioA = sin((1 - t) * halfTheta) / sinHalfTheta;
		float  ratioB = sin(t * halfTheta) / sinHalfTheta;
		//calculate Quaternion.
		res.w() = (qa.w() * ratioA + qb.w() * ratioB);
		res.x() = (qa.x() * ratioA + qb.x() * ratioB);
		res.y() = (qa.y() * ratioA + qb.y() * ratioB);
		res.z() = (qa.z() * ratioA + qb.z() * ratioB);
		return res;
	}
};
