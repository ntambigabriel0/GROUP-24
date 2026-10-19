### 1 How it was built?
A base class was established to hold shared atttributes and virtual functions, allowing derived classes to 
inherit common behaviour and override specific operations.
Private acccess specifiers and public accessor methods were implemented to secure internal data structures.
Finally, virtual functions and runtime dispatch were utilised to handle distinct object behaviours dynamically
through base class pointers and references.
The logic was validated by instatiating base and derived objects, confirming proper constructor / destructor chains, verifying data hiding boundaries and testing correct dynamic method binding.
There is built in safety through const correctness, virtual destructors to prevent memory leaks during polymorphic cleanup and process access encapsulation.

### 2 What classes were created?
Two primary classes were created: a base class e.g Shape or Account defining the common interface and protected attributes and one or more derived classes e.g Circle or SavingsAccount inheriting from the base class to implement specialised behaviours and polymorphic overrides.

### 3 What data does it hold?
The classes hold member variables categorised by access levels. Private variables such as raw coordinate buffers, balance amounts or identification numbers store core internal state while protected variables store shared attributes accessible to derived classes.

### 4 What functions operate on that data?
Public member functions operate on the class data. Additionally, virtual functions declared in the base class and overriddenin derived classes execute polymorphic operations based on the runtime object type.

### 5 Why some members are public and others hidden?
Member variables and internal helper logic are hidden to prevent unintended state corruption, enforce class invariants and maintain strict encapsulation. Only the necessary interface functions are kept public to expose safe interaction points to the  rest of the program.

### 6 How the rest of the program uses it?
BaseClass* obj = new DerivedClass();
obj->performAction(); // Resolves polymorphically at runtime
delete obj;

It integrates cleanly with standard C++ memory management and container libraries. Itcan be embedded into larger enterprise applications, game engienes, modular software systems requiring extensible architectures.

### 7 What issues were faced?
Object slicing when storing derived objects directly inside base class containers by value and memory leaks caused by missing virtual destructors in polymorphic base classes.
Future improvements include utilising smart pointers for safe resource management and expanding abstract interface hierachies.