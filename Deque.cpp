#include<iostream>
#include<vector>


template<typename T>
class Deque {

private:
	size_t arr_capacity;
	size_t bucket_size = 8;
	size_t deq_size = 0;
	T** arr;

public:
	template <bool is_const>
	struct common_iterator {
		T** arr;
		size_t capacity;
		int x;
		int y;
		T* ptr;

		common_iterator(void) = default;

		common_iterator(const common_iterator<true>& it) : arr(it.arr), capacity(it.capacity), x(it.x), y(it.y), ptr(it.ptr) {}
    
		common_iterator(const common_iterator<false>& it) : arr(it.arr), capacity(it.capacity), x(it.x), y(it.y), ptr(it.ptr) {}
    
		common_iterator(T** _arr, size_t _capacity, int _x, int _y) : arr(_arr), capacity(_capacity), x(_x), y(_y) {
			ptr = arr[x] + y;
		}

		common_iterator& operator=(const common_iterator& it) {
			arr = it.arr;
			capacity = it.capacity;
			x = it.x;
			y = it.y;
			ptr = it.ptr;
			return *this;
		}

		common_iterator& operator++() {
			++y;
			if (y == 8) {
				++x;
				y = 0;
			}
			ptr = arr[x] + y;
			return *this;
		}
    
		common_iterator& operator--() {
			--y;
			if (y == -1) {
				--x;
				y = 7;
			}
			ptr = arr[x] + y;
			return *this;
		}
    
		common_iterator& operator+=(int step) {
			step += y;
			int jump = step / 8;
			int ost = step % 8;
			x += jump;
			y = ost;
			ptr = arr[x] + y;
			return *this;
		}
    
		common_iterator& operator-=(int step) {
			step += (8 - y);
			int jump = (step - 1) / 8;
			int ost = (8 - (step % 8)) % 8;
			x -= jump;
			y = ost;
			ptr = arr[x] + y;
			return *this;
		}
    
		common_iterator operator+(int step) {
			common_iterator copy = *this;
			copy += step;
			return copy;
		}
    
		common_iterator operator-(int step) {
			common_iterator copy = *this;
			copy -= step;
			return copy;
		}

		bool operator<(const common_iterator& it) {
			return (x < it.x) || ((x == it.x) && (y < it.y));
		}
    
		bool operator>(const common_iterator& it) {
			return (x > it.x) || ((x == it.x) && (y > it.y));
		}
    
		bool operator<=(const common_iterator& it) {
			return (x < it.x) || ((x == it.x) && (y <= it.y));
		}
    
		bool operator>=(const common_iterator& it) {
			return (x > it.x) || ((x == it.x) && (y >= it.y));
		}
    
		bool operator==(const common_iterator& it) const{
			return (ptr == it.ptr);
		}
    
		bool operator!=(const common_iterator& it) const{
			return !(ptr == it.ptr);
		}
    
		int operator-(const common_iterator& it) const {
			return (x - it.x) * 8 + (y - it.y);
		}
    
		std::conditional_t<is_const, const T&, T&> operator*() {
			return *ptr;
		}
    
		std::conditional_t<is_const, const T*, T*> operator->() {
			return ptr;
		}

	};
	
	using iterator = common_iterator<false>;
	using const_iterator = common_iterator<true>;

	using reverse_iterator = std::reverse_iterator<iterator>;
	using const_reverse_iterator = std::reverse_iterator<const_iterator>;


	iterator _begin;
	iterator _end;


	Deque(void) : arr_capacity(1), arr(new T* [arr_capacity]) {
		arr[0] = reinterpret_cast<T*>(new int8_t[sizeof(T) * bucket_size]);
		_begin = iterator(arr, arr_capacity, 0, 0);
		_end = iterator(arr, arr_capacity, 0, 0);
	}
  
	Deque(const Deque& d) : arr_capacity(d.arr_capacity), deq_size(d.deq_size) {
		arr = new T * [arr_capacity];
		for (size_t i = 0; i != arr_capacity; i++)
			arr[i] = reinterpret_cast<T*>(new int8_t[sizeof(T) * bucket_size]);
		_begin = iterator(arr, arr_capacity, d.begin().x, d.begin().y);
		_end = iterator(arr, arr_capacity, d.end().x, d.end().y);
		size_t i = 0;
		try {
			for (; i != deq_size; ++i) {
				T value = *(d.begin() + i);
				new((_begin + i).ptr) T(value);
			}
		}
		catch (...) {
			for (size_t j = 0; j != i; ++j) {
				(_begin + i).ptr->~T();
			}
			for (size_t i = 0; i != arr_capacity; ++i)
				delete[] reinterpret_cast<int8_t*>(arr[i]);
			delete[] arr;
		}
	}
  
	Deque(int a) : Deque(a, T()) {}

	Deque(int a, const T& value) : arr_capacity((size_t)((a / 8) + 1) * 3), deq_size((size_t)(a)) {
		arr = new T * [arr_capacity];
		for (size_t i = 0; i != arr_capacity; i++)
			arr[i] = reinterpret_cast<T*>(new int8_t[sizeof(T) * bucket_size]);

		size_t i = 0;
		_begin = iterator(arr, arr_capacity, arr_capacity / 3, 0);
		_end = _begin + deq_size;

		try {
			for (; i != deq_size; ++i) {
				iterator s = _begin + i;
				new(s.ptr) T(value);
			}
		}
		catch (...) {
			for (size_t j = 0; j != i; ++j) {
				(_begin + i).ptr->~T();
			}
			for (size_t i = 0; i != arr_capacity; ++i)
				delete[] reinterpret_cast<int8_t*>(arr[i]);
			delete[] arr;
		}
	}

	Deque& operator=(const Deque& d) {
		for (size_t i = 0; i != deq_size; ++i) {
			(_begin + i).ptr->~T();
		}
		for (size_t i = 0; i != arr_capacity; ++i)
			delete[] reinterpret_cast<int8_t*>(arr[i]);
		delete[] arr;
		arr_capacity = d.arr_capacity;
		deq_size = d.deq_size;

		arr = new T * [arr_capacity];
		for (size_t i = 0; i != arr_capacity; i++)
			arr[i] = reinterpret_cast<T*>(new int8_t[sizeof(T) * bucket_size]);

		size_t i = 0;
		_begin = iterator(arr, arr_capacity, d.begin().x, d.begin().y);
		_end = _begin + deq_size;

		try {
			for (; i != deq_size; ++i) {
				new((_begin + i).ptr) T(*(d.begin() + i));
			}
		}
		catch (...) {
			for (size_t j = 0; j != i; ++j) {
				(_begin + i).ptr->~T();
			}
			for (size_t i = 0; i != arr_capacity; ++i)
				delete[] reinterpret_cast<int8_t*>(arr[i]);
			delete[] arr;
		}
		return *this;
	}

	size_t size() const {
		return deq_size;
	}

	T& operator[](size_t index) {
		iterator it = _begin;
		it += index;
		return *it;
	}
	T& at(size_t index) {
		iterator it = _begin;
		it += index;
		if (index >= deq_size) throw std::out_of_range("...");
		return *it;
	}

	const T& operator[](size_t index) const {
		iterator it = _begin;
		it += index;
		return *it;
	}
	const T& at(size_t index) const {
		iterator it = _begin;
		it += index;
		if (index >= deq_size) throw std::out_of_range("...");
		return *it;
	}

	void reserve() {
		std::cerr << "reserve" << ' ';
		size_t new_capacity = arr_capacity * 3;

		T** new_arr = new T * [new_capacity];
		for (size_t i = 0; i != new_capacity; i++)
			new_arr[i] = reinterpret_cast<T*>(new int8_t[sizeof(T) * bucket_size]);

		size_t i = 0;
		iterator new_begin = iterator(new_arr, new_capacity, new_capacity / 3, 0);

		try {
			for (; i != deq_size; ++i) {
				new((new_begin + i).ptr) T(*(_begin + i));
			}
		}
		catch (...) {
			for (size_t j = 0; j != i; ++j) {
				(new_begin + i).ptr->~T();
			}
			for (size_t i = 0; i != new_capacity; ++i)
				delete[] reinterpret_cast<int8_t*>(new_arr[i]);
			delete[] new_arr;
		}
		delete[] arr;
		arr = new_arr;
		arr_capacity = new_capacity;
		_begin = new_begin;
		_end = _begin + deq_size;
	}

	void push_back(const T& value) {
		if ((size_t)_end.x == arr_capacity) reserve();
		new(_end.ptr) T(value);
		++_end;
		++deq_size;
	}
  
	void pop_back() {
		(_end - 1).ptr->~T();
		--_end;
		--deq_size;
	}
  
	void push_front(const T& value) {
		if (_begin.x == 0 && _begin.y == 0) reserve();
		--_begin;
		new(_begin.ptr) T(value);
		++deq_size;
	}
  
	void pop_front() {
		_begin.ptr->~T();
		++_begin;
		--deq_size;
	}

	iterator begin() {
		return _begin;
	}
  
	iterator end() {
		return _end;
	}
  
	const_iterator cbegin() const {
		return const_iterator(_begin);
	}
  
	const_iterator cend() const {
		return const_iterator(_end);
	}
  
	const_iterator begin() const {
		return const_iterator(_begin);
	}
  
	const_iterator end() const {
		return const_iterator(_end);
	}
  
	reverse_iterator rbegin() {
		return reverse_iterator(_end - 1);
	}
  
	reverse_iterator rend() {
		return reverse_iterator(_begin() - 1);
	}
  
	const_reverse_iterator crbegin() const {
		return const_reverse_iterator(_end - 1);
	}
  
	const_reverse_iterator crend() const {
		return const_reverse_iterator(_begin - 1);
	}
  
	const_reverse_iterator rbegin() const {
		return const_reverse_iterator(_end - 1);
	}
  
	const_reverse_iterator rend() const {
		return const_reverse_iterator(_begin - 1);
	}

	void insert(iterator it, const T& value) {
		iterator iter = _end - 1;
		push_back(*iter);
		while (iter > it) {
			auto it_1 = iter - 1;
			*(iter) = *(it_1);
			--iter;
		}
		*(iter) = value;
	}
	void erase(iterator it) {
		iterator iter = it;
		while (iter < _end - 1) {
			*(iter) = *(iter + 1);
			++iter;
		}
		pop_back();
	}

	~Deque(void) {
		for (size_t i = 0; i != deq_size; ++i) {
			(_begin + i).ptr->~T();
		}
		for (size_t i = 0; i != arr_capacity; ++i)
			delete[] reinterpret_cast<int8_t*>(arr[i]);
		delete[] arr;
	}
  
};
