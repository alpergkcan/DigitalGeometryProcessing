#include "Data.h"
// VEC3
	Vec3::Vec3() : x(0), y(0), z(0) {};
        Vec3::Vec3(const float n) : x(n), y(n), z(n) {};
	Vec3::Vec3(const float x1, const float y1, const float z1) : x(x1), y(y1), z(z1) {};

	Vec3 Vec3::normalized() {
		float mag = sqrt(x * x + y * y + z * z);
		return Vec3(x / mag, y / mag, z / mag);
	}
	void Vec3::normalize() {
		float mag = sqrt(x * x + y * y + z * z);

		x /= mag;
		y /= mag;
		z /= mag;
	}

	float Vec3::magnitude() {
		return sqrt(x * x + y * y + z * z);
	}

	float* Vec3::toFloat3() {
		float* array = new float[3];

		array[0] = x;
		array[1] = y;
		array[2] = z;

		return array;
	}
	std::ostream& operator<<(std::ostream& os, Vec3 rhs) {
		return os << "[ " << rhs.x << "  " << rhs.y << "  " << rhs.z << " ]";
	}

	// OPERATORS
	Vec3 Vec3::operator/(const float& n) const { // Scaler Div 
		return Vec3(x / n, y / n, z / n);
	}
	Vec3 Vec3::operator/(const int& n) const { // Scaler Div 
		return Vec3(x / n, y / n, z / n);
	}
	Vec3 Vec3::operator/(const uint& n) const { // Scaler Div 
		return Vec3(x / n, y / n, z / n);
	}
	Vec3 operator/(const float& n, const Vec3& vec) {
		return Vec3(n / vec.x, n / vec.y, n / vec.z);
	}
	Vec3 Vec3::operator^(const float& n) const { // Scaler Power
		return Vec3(pow(x, n), pow(y, n), pow(z, n));
	}
	Vec3 operator^(const float& n, const Vec3& vec) { // Scaler Power
		return Vec3(pow(n, vec.x), pow(n, vec.y), pow(n, vec.z));
	}
	Vec3 Vec3::operator*(const float& n) const { // Scalar Prod
		return Vec3(x * n, y * n, z * n);
	}
	Vec3 operator*(const float& n, const Vec3& vec) { // Scalar Prod
		return Vec3(vec.x * n, vec.y * n, vec.z * n);
	}
	Vec3 Vec3::operator+(const float& n) const { // vec addition
		return Vec3(x + n, y + n, z + n);
	}
	Vec3 operator+(const float& n, const Vec3& vec) { // vec addition
		return Vec3(vec.x + n, vec.y + n, vec.z + n);
	}
	Vec3 Vec3::operator-(const float& n) const { // vec subs
		return Vec3(x - n, y - n, z - n);
	}
	Vec3 operator-(const float& n, const Vec3& vec) { // vec subs
		return Vec3(n - vec.x, n - vec.y, n - vec.z);
	}


	bool Vec3::operator==(const Vec3& rhs) const {
		return (x == rhs.x && y == rhs.y && z == rhs.z);
	}


	Vec3 Vec3::operator*(const Vec3 rhs) const { // Dot
		return Vec3(x * rhs.x, y * rhs.y, z * rhs.z);
	}
	Vec3 Vec3::operator/(const Vec3 rhs) const { // divide
		return Vec3(x / rhs.x, y / rhs.y, z / rhs.z);
	}
	Vec3 Vec3::operator-(const Vec3& rhs) const { // vec subs
		return Vec3(x - rhs.x, y - rhs.y, z - rhs.z);
	}
	Vec3 Vec3::operator+(const Vec3& rhs) const { // vec addition
		return Vec3(x + rhs.x, y + rhs.y, z + rhs.z);
	}


	float& Vec3::operator[](int idx) {
		switch (idx) {
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
		default:
			assert(idx < 0 or idx > 2);
			break;
		}
	}
	float Vec3::operator[](int idx) const {
		switch (idx) {
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
		default:
			assert(idx < 0 or idx > 2);
			return 0;
			break;
		}
	}
// END OF VEC3



// VEC 4
	Vec4::Vec4() : x(0), y(0), z(0), w(0) {};
	Vec4::Vec4(const float n) : x(n), y(n), z(n), w(n) {};
	Vec4::Vec4(Vec3 cpy, const float n) : x(cpy.x), y(cpy.y), z(cpy.z), w(n) {};
	Vec4::Vec4(const float n, Vec3 cpy) : y(cpy.x), z(cpy.y), w(cpy.z), x(n) {};
	Vec4::Vec4(const float x1, const float y1, const float z1, const float w1) : x(x1), y(y1), z(z1), w(w1) {};

	Vec4 Vec4::normalized() {
		float mag = sqrt(x * x + y * y + z * z + w * w);
		return Vec4(x / mag, y / mag, z / mag, w / mag);
	}
	void Vec4::normalize() {
		float mag = sqrt(x * x + y * y + z * z + w * w);

		x /= mag;
		y /= mag;
		z /= mag;
		w /= mag;
	}

	float Vec4::magnitude() {
		return sqrt(x * x + y * y + z * z);
	}

	float* Vec4::toFloat4() {
		float* array = new float[4];

		array[0] = x;
		array[1] = y;
		array[2] = z;
		array[3] = w;

		return array;
	}

	std::ostream& operator<<(std::ostream& os, Vec4 rhs) {
		return os << "[ " << rhs.x << "  " << rhs.y << "  " << rhs.z << "  " << rhs.w << " ]";
	}

	// OPERATORS
	Vec4 Vec4::operator/(const float& n) const { // Scaler Div 
		return Vec4(x / n, y / n, z / n, w / n);
	}
	Vec4 operator/(const float& n, const Vec4& vec) {
		return Vec4(n / vec.x, n / vec.y, n / vec.z, n / vec.w);
	}
	Vec4 Vec4::operator^(const float& n) const { // Scaler Power
		return Vec4(pow(x, n), pow(y, n), pow(z, n), pow(w, n));
	}
	Vec4 operator^(const float& n, const  Vec4& vec) { // Scaler Power
		return Vec4(pow(n, vec.x), pow(n, vec.y), pow(n, vec.z), pow(n, vec.w));
	}
	Vec4 Vec4::operator*(const float& n) const { // Scalar Prod
		return Vec4(x * n, y * n, z * n, w * n);
	}
	Vec4 operator*(const float& n, const  Vec4& vec) { // Scalar Prod
		return Vec4(vec.x * n, vec.y * n, vec.z * n, vec.w * n);
	}
	Vec4 Vec4::operator+(const float& n) const { // vec addition
		return Vec4(x + n, y + n, z + n, w + n);
	}
	Vec4 operator+(const float& n, const Vec4& vec) { // vec addition
		return Vec4(vec.x + n, vec.y + n, vec.z + n, vec.w + n);
	}
	Vec4 Vec4::operator-(const float& n) const { // vec subs
		return Vec4(x - n, y - n, z - n, w - n);
	}
	Vec4 operator-(const float& n, const Vec4& vec) { // vec subs
		return Vec4(n - vec.x, n - vec.y, n - vec.z, n - vec.z);
	}


	Vec4 Vec4::operator*(const Vec4& rhs) const { // Dot
		return Vec4(x * rhs.x, y * rhs.y, z * rhs.z, w * rhs.w);
	}
	Vec4 Vec4::operator/(const Vec4& rhs) const { // divide
		return Vec4(x / rhs.x, y / rhs.y, z / rhs.z, w / rhs.w);
	}
	Vec4 Vec4::operator-(const Vec4& rhs) const { // vec subs
		return Vec4(x - rhs.x, y - rhs.y, z - rhs.z, w - rhs.w);
	}
	Vec4 Vec4::operator+(const Vec4& rhs) const { // vec addition
		return Vec4(x + rhs.x, y + rhs.y, z + rhs.z, w + rhs.z);
	}

	float& Vec4::operator[](int idx) {
		switch (idx) {
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
		case 3:
			return w;
		default:
			assert(idx < 0 or idx > 3);
			break;
		}
	}
	float Vec4::operator[](int idx) const {
		switch (idx) {
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
		case 3:
			return w;
		default:
			assert(idx < 0 or idx > 3);
			return 0;
			break;
		}
	}
// END OF VEC 4


// MAT 3
        Mat3::Mat3() : x(Vec3(1,0,0)), y(Vec3(0,1,0)), z(Vec3(0,0,1)) {}
	Mat3::Mat3(const float n) : x(Vec3(n)), y(Vec3(n)), z(Vec3(n)) {}

	Mat3::Mat3(const Vec3& _x, const Vec3& _y, const Vec3& _z) : x(_x), y(_y), z(_z) {}

	Mat3::Mat3(const float*& array) : x(Vec3(array[0], array[1], array[2])),
		y(Vec3(array[3], array[4], array[5])),
		z(Vec3(array[6], array[7], array[8])) {}

	Mat3::Mat3( float**& array) :      x(Vec3(array[0][0], array[0][1], array[0][2])),
									   y(Vec3(array[1][0], array[1][1], array[1][2])),
									   z(Vec3(array[2][0], array[2][1], array[2][2])) {}

	Mat3::Mat3(vector<float>& array) : x(Vec3(array[0], array[1], array[2])),
		y(Vec3(array[3], array[4], array[5])),
		z(Vec3(array[6], array[7], array[8])) {}

	std::ostream& operator<<(std::ostream& os, const Mat3& mat) {
		os << std::setprecision(4) << "[ ";
		for (int i = 0; i < 2; i++)
			os << "  " << mat[i] << "\n";
		return os << "  " << mat[2] << " ]";
	}

	Mat3 Mat3::operator+(const float& n) const {
		return Mat3(x + n, y + n, z + n);
	}
	Mat3 operator+(const float& n, const Mat3& mat) {
		return Mat3(n + mat.x, n + mat.y, n + mat.z);
	}
	Mat3 Mat3::operator-(const float& n) const {
		return Mat3(x - n, y - n, z - n);
	}
	Mat3 operator-(const float& n, const Mat3& mat) {
		return Mat3(n - mat.x, n - mat.y, n - mat.z);
	}
	Mat3 Mat3::operator*(const float& n) const {
		return Mat3(x * n, y * n, z * n);
	}
	Mat3 operator*(const float& n, const Mat3& mat) {
		return Mat3(n * mat.x, n * mat.y, n * mat.z);
	}
	Mat3 Mat3::operator/(const float& n) const {
		return Mat3(x / n, y / n, z / n);
	}

	Mat3 operator/(const float& n, const Mat3& mat) {
		return Mat3(n / mat.x, n / mat.y, n / mat.z);
	}




	Mat3 Mat3::operator+(const Mat3& rhs) const {
		return Mat3(x + rhs.x, y + rhs.y, z + rhs.z);
	}
	Mat3 Mat3::operator-(const Mat3& rhs) const {
		return Mat3(x - rhs.x, y - rhs.y, z - rhs.z);
	}
	Mat3 Mat3::operator*(const Mat3& rhs) const {
		return Mat3(x * rhs.x, y * rhs.y, z * rhs.z);
	}
	Vec3 Mat3::operator*(const Vec3& vec) const {
		return Vec3(Data::dot(x, vec), Data::dot(y, vec), Data::dot(z, vec));
	}
	Vec3 operator*(const Vec3& vec, const Mat3& _mat) {
		Mat3 mat = Data::transpose(_mat);
		return Vec3(Data::dot(vec, mat.x), Data::dot(vec, mat.y), Data::dot(vec, mat.z));
	}
	Mat3 Mat3::operator/(const Mat3& rhs) const {
		return Mat3(x / rhs.x, y / rhs.y, z / rhs.z);
	}

	Vec3& Mat3::operator[](int idx) {
		switch (idx) {
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
		default:
			assert(idx < 0 or idx > 2);
			break;
		}
	}
	Vec3 Mat3::operator[](int idx) const {
		switch (idx) {
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
		default:
			assert(idx < 0 or idx > 2);
			return 0;
			break;
		}
	}
	float* Mat3::toFloat3() {
		float* res = new float[12];

		for (size_t i = 0; i < 3; i++)
		{
			res[i] = x[i];
			res[i + 4] = y[i];
			res[i + 8] = z[i];
		}
		return res;
	}
	float** Mat3::toFloat4_2() {
		float** res = new float* [4];

		for (size_t i = 0; i < 4; i++)
		{
			res[i] = new float[4];

			res[i][0] = x.x;
			res[i][1] = x.y;
			res[i][2] = x.z;
			res[i][3] = 0;
		}
		return res;
	}

// END OF MAT 3



// MAT 4
	Mat4::Mat4() : x(Vec4()), y(Vec4()), z(Vec4()), w(Vec4()) {}
	Mat4::Mat4(const float n) : x(Vec4(n)), y(Vec4(n)), z(Vec4(n)), w(Vec4(n)) {}

	Mat4::Mat4(const Vec4& _x, const Vec4& _y, const Vec4& _z, const Vec4& _w) : x(_x), y(_y), z(_z), w(_w) {}


	Mat4::Mat4(Mat3& x3, Vec3& col, Vec3& row, const float deep) : x(x3.x.x, x3.x.y, x3.x.z, col.x),
		y(x3.y.x, x3.y.y, x3.y.z, col.y),
		z(x3.z.x, x3.z.y, x3.z.z, col.z),
		w(row.x, row.y, row.z, deep) {}

	Mat4::Mat4( float* array) : x(Vec4(array[0], array[1], array[2], array[3])),
		y(Vec4(array[4], array[5], array[6], array[7])),
		z(Vec4(array[8], array[9], array[10], array[11])),
		w(Vec4(array[12], array[13], array[14], array[15])) {}
	
	Mat4::Mat4( float** array) :       x(Vec4(array[0][0], array[0][1], array[0][2], array[0][3])),
								       y(Vec4(array[1][0], array[1][1], array[1][2], array[1][3])),
									   z(Vec4(array[2][0], array[2][1], array[2][2], array[2][3])),
									   w(Vec4(array[3][0], array[3][1], array[3][2], array[3][3]))	{}

	Mat4::Mat4(vector<float>& array) :  x(Vec4(array[0], array[1], array[2], array[3])),
										y(Vec4(array[4], array[5], array[6], array[7])),
										z(Vec4(array[8], array[9], array[10], array[11])),
										w(Vec4(array[12], array[13], array[14], array[15])) {}

	std::ostream& operator<<(std::ostream& os, const Mat4& mat) {
		os << std::setprecision(4) << "[ ";
		for (int i = 0; i < 3; i++)
			os << "  " << mat[i] << "\n";
		return os << "  " << mat[3] << " ]";
	}
	Mat4 Mat4::operator+(const float& n) const {
		return Mat4(x + n, y + n, z + n, w + n);
	}
	Mat4 operator+(const float& n, const Mat4& mat) {
		return Mat4(mat.x + n, mat.y + n, mat.z + n, mat.w + n);
	}
	Mat4 Mat4::operator-(const float& n) const {
		return Mat4(x - n, y - n, z - n, w - n);
	}
	Mat4 operator-(const float& n, const Mat4& mat) {
		return Mat4(n - mat.x, n - mat.y, n - mat.z, n - mat.w);
	}
	Mat4 Mat4::operator*(const float& n) const {
		return Mat4(x * n, y * n, z * n, w * n);
	}
	Mat4 operator*(const float& n, const Mat4& mat) {
		return Mat4(mat.x * n, mat.y * n, mat.z * n, mat.w * n);
	}
	Mat4 Mat4::operator/(const float& n) const {
		return Mat4(x / n, y / n, z / n, w / n);
	}
	Mat4 operator/(const float& n, const Mat4& mat) {
		return Mat4(n / mat.x, n / mat.y, n / mat.z, n / mat.w);
	}



	Mat4 Mat4::operator+(const Mat4& rhs) const {
		return Mat4(x + rhs.x, y + rhs.y, z + rhs.z, w + rhs.w);
	}
	Mat4 Mat4::operator-(const Mat4& rhs) const {
		return Mat4(x - rhs.x, y - rhs.y, z - rhs.z, w - rhs.w);
	}
	Mat4 Mat4::operator*(const Mat4& rhs) const {
		return Mat4(x * rhs.x, y * rhs.y, z * rhs.z, w * rhs.w);
	}
	Vec4 Mat4::operator*(const Vec4& vec) const {
		return Vec4(Data::dot(x, vec), Data::dot(y, vec), Data::dot(z, vec), Data::dot(w, vec));
	}
	Vec4 operator*(const Vec4& vec, const Mat4& _mat) {
		Mat4 mat = Data::transpose(_mat);
		return Vec4(Data::dot(vec, mat.x), Data::dot(vec, mat.y), Data::dot(vec, mat.z), Data::dot(vec, mat.w));
	}
	Mat4 Mat4::operator/(const Mat4& rhs) const {
		return Mat4(x / rhs.x, y / rhs.y, z / rhs.z, w / rhs.w);
	}




	Vec4& Mat4::operator[](int idx) {
		switch (idx) {
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
		case 3:
			return w;
		default:
			assert(idx < 0 or idx > 3);
			break;
		}
	}
	Vec4 Mat4::operator[](int idx) const {
		switch (idx) {
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
		case 3:
			return w;
		default:
			assert(idx < 0 or idx > 3);
			return 0;
			break;
		}
	}

	float* Mat4::toFloat4() {
		float* res = new float[16];

		for (size_t i = 0; i < 4; i++)
		{
			res[i]    = x[i];
			res[i+4]  = y[i];
			res[i+8]  = z[i];
			res[i+12] = w[i];
		}
		return res;
	}
	float** Mat4::toFloat4_2() {
		float** res = new float*[4];

		for (size_t i = 0; i < 4; i++)
		{
			res[i] = new float[4];

			res[i][0] = x.x;
			res[i][1] = x.y;
			res[i][2] = x.z;
			res[i][3] = x.w;
		}
		return res;
	}

// MAT 4


// DATA
	Mat3 Data::outer(const Vec3& c, const Vec3& r) {
		return Mat3(c * r.x, c * r.y, c * r.z);
	}
	Mat4 Data::outer(const Vec4& c, const Vec4& r) {
		return Mat4(c * r.x, c * r.y, c * r.z, c * r.w);
	}


	Vec3 Data::cross(const Vec3& lhs, const Vec3& rhs) { // Cross
		return Vec3((lhs.y * rhs.z) - (lhs.z * rhs.y), (lhs.z * rhs.x) - (lhs.x * rhs.z), (lhs.x * rhs.y) - (lhs.y * rhs.x));
	}


	float Data::dot(const Vec3& lhs, const Vec3& rhs) { // Dot
		return ((lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z));
	}
	float Data::dot(const Vec4& lhs, const Vec4& rhs) { // Dot
		return ((lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z) + (lhs.w * rhs.w));
	}


	Mat4 Data::transpose(const Mat4& mat) {
		return Mat4(Vec4(mat.x.x, mat.y.x, mat.z.x, mat.w.x),
			Vec4(mat.x.y, mat.y.y, mat.z.y, mat.w.y),
			Vec4(mat.x.z, mat.y.z, mat.z.z, mat.w.z),
			Vec4(mat.x.w, mat.y.w, mat.z.w, mat.w.w));
	}
	Mat3 Data::transpose(const Mat3& mat) {
		return Mat3(Vec3(mat.x.x, mat.y.x, mat.z.x),
			Vec3(mat.x.y, mat.y.y, mat.z.y),
			Vec3(mat.x.z, mat.y.z, mat.z.z));
	}

	template<class T>
	T Data::getCofactor(T A, int p, int q, int n)
	{
		int i = 0, j = 0;
		T temp;
		// Looping for each element of the matrix
		for (int row = 0; row < n; row++)
		{
			for (int col = 0; col < n; col++)
			{
				//  Copying into temporary matrix only those element
				//  which are not in given row and column
				if (row != p && col != q)
				{
					temp[i][j++] = A[row][col];

					// Row is filled, so increase row index and
					// reset col index
					if (j == n - 1)
					{
						j = 0;
						i++;
					}
				}
			}
		}
		return temp;
	}

	/* Recursive function for finding determinant of matrix.
		n is current dimension of A[][]. */
	template <class T>
	float Data::determinant(T mat, int n, int size)
	{
		//  Base case : if matrix contains single element
		if (n == 1)
			return mat[0][0];

		int D = 0; // Initialize result

		T temp;
		int sign = 1;  // To store sign multiplier

			// Iterate for each element of first row
		for (int f = 0; f < n; f++)
		{
			// Getting Cofactor of A[0][f]
			temp = getCofactor(mat, 0, f, n); // CHANGE!!
			D += sign * mat[0][f] * determinant(temp, n - 1, size);

			// terms are to be added with alternate sign
			sign = -sign;
		}
		return D;
	}

	// Function to get adjoint of A[N][N] in adj[N][N].
	template <class T>
	T Data::adjoint(T mat, int size)
	{
		// temp is used to store cofactors of A[][]
		int sign = 1;
		T temp;
		T adj;

		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				// Get cofactor of A[i][j]
				temp = getCofactor(mat, i, j, size); // CHANGE!!

				// sign of adj[j][i] positive if sum of row
				// and column indexes is even.
				sign = ((i + j) % 2 == 0) ? 1 : -1;

				// Interchanging rows and columns to get the
				// transpose of the cofactor matrix
				adj[j][i] = (sign) * (determinant(temp, size - 1, size));
			}
		}
		return adj;
	}



	mat3 Data::inverse( mat3& mat)
	{
		// Find determinant of A[][]
		float det = determinant(mat, 3, 3); // CHANGE!
		if (det == 0)
		{
			std::cout << "Singular matrix, can't find its inverse";
			assert(0);
		}

		// Find adjoint
		mat3 adj = adjoint(mat, 3);
		mat3 inverse;

		// Find Inverse using formula "inverse(A) = adj(A)/det(A)"
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				inverse[i][j] = adj[i][j] / (float)det;

		return inverse;
	}
	mat4 Data::inverse( mat4& mat)
	{
		// Find determinant of A[][]
		float det = determinant(mat, 4, 4); // CHANGE!
		if (det == 0)
		{
			std::cout << "Singular matrix, can't find its inverse";
			assert(0);
		}

		// Find adjoint
		mat4 adj = adjoint(mat, 4);
		mat4 inverse;

		// Find Inverse using formula "inverse(A) = adj(A)/det(A)"
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				inverse[i][j] = adj[i][j] / (float)det;

		return inverse;
	}
	
// DATA
	
	float Quaternion::magnitude() {
		return pow(w, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2);
	}
	void Quaternion::normalize() {
		float mag = magnitude();
		w /= mag;
		x /= mag;
		y /= mag;
		z /= mag;
	}
	Quaternion Quaternion::normalized(Quaternion qua) {
		Quaternion res;
		float mag = qua.magnitude();
		res.w = qua.w /mag;
		res.x = qua.x /mag;
		res.y = qua.y /mag;
		res.z = qua.z /mag;
		return res;
	}

	Mat3 Quaternion::to_mat(Quaternion& quat) 
	{
		if (!(1 - FLT_EPSILON <= quat.magnitude() <= 1 + FLT_EPSILON)) {

		}


		float p,r,s,k, l, m;
		


		p = 2 * quat.x * quat.y;
		r = 2 * quat.x * quat.z;
		s = 2 * quat.y * quat.z;


		k = 2 * quat.w * quat.z;
		l = 2 * quat.w * quat.y;
		m = 2 * quat.w * quat.x;

		return Mat3(vec3(1 - 2 * pow(quat.y, 2) - 2 * pow(quat.z, 2), p + k, r - l),
					vec3(p - k, 1 - 2 * pow(quat.x, 2) - 2 * pow(quat.z, 2), s + m),
					vec3(r + l, s - m, 1 - 2 * pow(quat.x, 2) - 2 * pow(quat.y, 2)) );
	}

	Quaternion Quaternion::to_quat(Mat3& mat)
	{
		float sq = 0.25f * (1 + mat[1][1] + mat[2][2] + mat[3][3]);
		float w, x, y, z;
		if (sq > FLT_EPSILON)
		{
	
			w = sqrt(sq);
			x = (mat[2][3] - mat[3][2]) / (4 * w);
			y = (mat[3][1] - mat[1][3]) / (4 * w);
			z = (mat[1][2] - mat[2][1]) / (4 * w);
		}else{
			w = 0;
			sq = -0.5f * (mat[2][2] + mat[3][3]);
			if (sq > FLT_EPSILON)
			{
				x = sqrt(sq);
				y = mat[1][2] / (2 * x);
				z = mat[1][3] / (2 * x);
			}else{
				x = 0;
				sq = 0.5f * (1 - mat[3][3]);
				if (sq > FLT_EPSILON)
				{
					y = sqrt(sq);
					z = mat[2][3] / (2 * y);
				}else{
					y = 0;
					z = 1;
				}
			}
		}
		return Quaternion(w,x,y,z);
	}
	
