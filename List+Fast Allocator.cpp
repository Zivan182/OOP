#include<iostream>


template<typename T, typename Alloc = std::allocator<T>>
class List {

private:

	struct Node {

		T value;
		Node* prev = nullptr;
		Node* next = nullptr;

		Node(const T& value, Node* prev, Node* next) : value(value), prev(prev), next(next) {}

		Node(Node* prev, Node* next) : prev(prev), next(next) {}
	};

	Node* fake = reinterpret_cast<Node*>(new int8_t[sizeof(Node)]);
	size_t sz = 0;

	using newAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<Node>;
	using AllocTraits = std::allocator_traits<newAlloc>;

	newAlloc alloc;

	template<bool IsConst>
	struct common_iterator {

		Node* ptr;

		common_iterator(Node* ptr) : ptr(ptr) {}

		common_iterator(const common_iterator<false>& iter) : ptr(iter.ptr) {}

		std::conditional_t<IsConst, const T&, T&> operator*() {
			return ptr->value;
		}

		common_iterator& operator++() {
			ptr = ptr->next;
			return *this;
		}

		common_iterator operator++(int) {
			ptr = ptr->next;
			return ptr->prev;
		}

		common_iterator& operator--() {
			ptr = ptr->prev;
			return *this;
		}

		common_iterator operator--(int) {
			ptr = ptr->prev;
			return ptr->next;
		}

		bool operator==(const common_iterator& iter) {
			return ptr == iter.ptr;
		}

		bool operator!=(const common_iterator& iter) {
			return ptr != iter.ptr;
		}

		using iterator_category = std::bidirectional_iterator_tag;
		using value_type = T;
		using pointer = std::conditional_t<IsConst, const T*, T*>;
		using reference = std::conditional_t<IsConst, const T&, T&>;
		using difference_type = std::ptrdiff_t;

	};

public:

	using iterator = common_iterator<false>;
	using const_iterator = common_iterator<true>;
	using reverse_iterator = std::reverse_iterator<iterator>;
	using const_reverse_iterator = std::reverse_iterator<const_iterator>;

	List() {
		fake->next = fake;
		fake->prev = fake;
	}

	explicit List(const Alloc& alloc) : alloc(alloc) {
		fake->next = fake;
		fake->prev = fake;
	}

	List(size_t count, const T& value, const Alloc& alloc = Alloc()) :
		alloc(AllocTraits::select_on_container_copy_construction(alloc)) {
		fake->next = fake;
		fake->prev = fake;
		for (size_t i = 0; i < count; ++i) {
			push_back(value);
		}
	}

	List(size_t count, const Alloc& _alloc = Alloc()) :
		alloc(AllocTraits::select_on_container_copy_construction(_alloc)) {
		fake->next = fake;
		fake->prev = fake;
		for (size_t i = 0; i < count; ++i) {
			push_back();
		}
	}

	List(const List& another) :
		alloc(AllocTraits::select_on_container_copy_construction(another.alloc)) {
		fake->next = fake;
		fake->prev = fake;
		iterator it = another.begin();
		for (size_t i = 0; i < another.sz; ++i, ++it) {
			push_back(*it);
		}
	}

	List& operator=(const List& another) {
		if (this == &another) {
			return *this;
		}
		while (sz != 0) {
			pop_back();
		}
		fake->next = fake;
		fake->prev = fake;
		if (AllocTraits::propagate_on_container_copy_assignment::value && alloc != another.alloc) {
			alloc = another.alloc;
		}
		iterator it = another.begin();
		for (size_t i = 0; i < another.sz; ++i, ++it) {
			push_back(*it);
		}
		return *this;
	}

	newAlloc get_allocator() const {
		return alloc;
	}

	iterator begin() const {
		return iterator(fake->next);
	}

	iterator end() const {
		return iterator(fake);
	}

	const_iterator cbegin() const {
		return const_iterator(fake->next);
	}

	const_iterator cend() const {
		return const_iterator(fake);
	}

	reverse_iterator rbegin() const {
		return reverse_iterator(fake);
	}

	reverse_iterator rend() const {
		return reverse_iterator(fake->prev);
	}

	const_reverse_iterator crbegin() const {
		return const_reverse_iterator(fake);
	}

	const_reverse_iterator crend() const {
		return const_reverse_iterator(fake->prev);
	}

	size_t size() const {
		return sz;
	}

	void push_back(const T& value) {
		Node* ptr = AllocTraits::allocate(alloc, 1);
		AllocTraits::construct(alloc, ptr, value, fake->prev, fake);
		(fake->prev)->next = ptr;
		fake->prev = ptr;
		++sz;
	}

	void push_back() {
		Node* ptr = AllocTraits::allocate(alloc, 1);
		AllocTraits::construct(alloc, ptr, fake->prev, fake);
		(fake->prev)->next = ptr;
		fake->prev = ptr;
		++sz;
	}

	void push_front(const T& value) {
		Node* ptr = AllocTraits::allocate(alloc, 1);
		AllocTraits::construct(alloc, ptr, value, fake, fake->next);
		(fake->next)->prev = ptr;
		fake->next = ptr;
		++sz;
	}

	void pop_back() {
		Node* ptr = fake->prev;
		((fake->prev)->prev)->next = fake;
		fake->prev = (fake->prev)->prev;
		AllocTraits::destroy(alloc, ptr);
		AllocTraits::deallocate(alloc, ptr, 1);
		--sz;
	}

	void pop_front() {
		Node* ptr = fake->next;
		((fake->next)->next)->prev = fake;
		fake->next = (fake->next)->next;
		AllocTraits::destroy(alloc, ptr);
		AllocTraits::deallocate(alloc, ptr, 1);
		--sz;
	}

	void insert(const_iterator iter, const T& value) {
		Node* ptr = AllocTraits::allocate(alloc, 1);
		AllocTraits::construct(alloc, ptr, value, (iter.ptr)->prev, iter.ptr);
		(iter.ptr->prev)->next = ptr;
		iter.ptr->prev = ptr;
		++sz;
	}

	void erase(const_iterator iter) {
		(iter.ptr->prev)->next = iter.ptr->next;
		(iter.ptr->next)->prev = iter.ptr->prev;
		AllocTraits::destroy(alloc, iter.ptr);
		AllocTraits::deallocate(alloc, iter.ptr, 1);
		--sz;
	}

	~List() {
		while (sz != 0) {
			pop_back();
		}
	}
};



template<size_t chunkSize>
class FixedAllocator {

private:

	static FixedAllocator* fixed;
	static const long long sz = 50000000;
	static int count;
	int8_t* pool;
	int number = 0;

	FixedAllocator() : pool(reinterpret_cast<int8_t*>(::operator new(sz))) {};

public:

	static FixedAllocator* GetAllocator() {
		++count;
		if (fixed) {
			return fixed;
		}
		fixed = new FixedAllocator();
		return fixed;
	}

	static void Destroy() {
		--count;
		if (count == 0) {
			fixed->number = 0;
		}
	}

	int8_t* allocate(size_t n) {
		number += n * chunkSize;
		return (pool + number - n * chunkSize);
	}

	void deallocate(void* ptr, size_t) {
		return;
	}

	~FixedAllocator() {
		::operator delete(pool);
	}
};

template<size_t chunkSize>
FixedAllocator<chunkSize>* FixedAllocator<chunkSize>::fixed = nullptr;

template<size_t chunkSize>
int FixedAllocator<chunkSize>::count = 0;



template<typename T>
class FastAllocator {

private:

	FixedAllocator<sizeof(T)>* fixed = FixedAllocator<sizeof(T)>::GetAllocator();

public:

	using value_type = T;

	FastAllocator(void) {}

	FastAllocator(const FastAllocator&) {}

	template<class U>
	FastAllocator(const FastAllocator<U>&) {}

	FastAllocator& operator=(const FastAllocator&) {
		return *this;
	}

	T* allocate(size_t n) {
		if (sizeof(T) <= 32) {
			return reinterpret_cast<T*>(fixed->allocate(n));
		}
		return reinterpret_cast<T*>(::operator new(n * sizeof(T)));
	}

	void deallocate(T* ptr, size_t) {
		if (sizeof(T) <= 32) {
			return;
		}
		::operator delete(ptr);
	}

	template<typename U>
	struct rebind {
		typedef FastAllocator<U> other;
	};

	~FastAllocator() {
		FixedAllocator<sizeof(T)>::Destroy();
	}
};

template<typename T>
bool operator==(const FastAllocator<T>&, const FastAllocator<T>&) {
	return true;
}

template<typename T>
bool operator!=(const FastAllocator<T>&, const FastAllocator<T>&) {
	return false;
}
