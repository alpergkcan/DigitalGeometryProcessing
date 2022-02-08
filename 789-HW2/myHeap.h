#pragma once

class distpair{
public:
	int x;
	float y;
	distpair(int x1 = -1, float y1 = 2e9): x(x1), y(y1){}
	
	bool operator<(distpair rhs) {
		return y < rhs.y;
	}
	
	bool operator==(distpair rhs) {
		return y == rhs.y;
	}
	
	bool operator>(distpair rhs) {
		return y > rhs.y;
	}
};


template <class T>
class minHeap {
private:
	T* array;
	int theSize;
	int capacity;
public:
	minHeap ( int cap = 100) {
		theSize = 0;
		capacity = cap;
		array = new T[capacity+1];
	}
	minHeap ( T* a = nullptr, int n=0) {
	 	theSize = n;
	 	capacity = n+1;
		
	 	array = new T[capacity];
		for (int i = 1; i <= n; ++i) 
			array[i] = a[i-1];
	 	buildHeap();
	}
		
	bool isEmpty( ) {
		return (theSize == 0);
	}
	
	T & getMin( ){
		return array[1];
	}

	//	void deleteMin( T & minItem );
	void makeEmpty( ){
		delete [] array;
	}

	void resize(int x ) {
		capacity = x;
		T* newArray = new T[capacity];
		
		for (int i = 0; i < theSize; ++i)
			newArray[i] = array[i];
		
		delete [] array;
		array = newArray;

	}
	
	void insert( T & x ) {
		array[ 0 ] = x;   // initialize sentinel
		if( theSize + 1 == capacity )
			resize( theSize * 2 + 1 );

		// Percolate up
		int hole = ++theSize;
		for( ; x < array[ hole / 2 ]; hole /= 2 )
			array[ hole ] = array[ hole / 2 ];
		array[ hole ] = x;
	}

	T deleteMin( ) {
		if( isEmpty( ) )
			throw "NON EMPTY HEAP";
		T temp = array[1];
		array[1] = array[ theSize-- ];
		percolateDown( 1 );
		return temp;
	}

	void percolateDown( int hole ) {
		int child;
		T tmp = array[ hole ];

		for( ; hole * 2 <= theSize; hole = child ) {
			child = hole * 2;
			
			if( child != theSize && array[child + 1] < array[child])
				child++;
			if( array[ child ] < tmp )
				array[ hole ] = array[ child ];
			else
				break;
		}
		
		array[ hole ] = tmp;
	}

	void buildHeap( ) {
		for( int i = theSize / 2; i > 0; i-- )
			percolateDown( i );
	}


	void findAndPercUp(int idx, float val){

		for (int i = 1; i <= theSize; ++i) {
			if(array[i].x==idx){
				array[i].y = val;
				
				// Percolate up

				distpair x = array[i];
				array[0] = x;
				for( ; x < array[ i / 2 ]; i /= 2 )
					array[ i ] = array[ i / 2 ];
				array[ i ] = x;
				return;
			}
				
		}

	}
};
// FIX!!!!
//
// BUILD HEAP
// DELETE MIN
// CONSTRUCTOR n
// isEMPTY()
