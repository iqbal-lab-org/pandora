# Testing guidelines

Here are described some conventions to write tests in `pandora`. While you are encouraged to follow these conventions,
we are not strict in the sense that you can bypass these conventions in case your test case needs, or if by following these
conventions your tests are hard to read. The most important thing is to make tests readable and easy to understand. Also,
we should consider tests as important as production code, because they are the ones enabling us to change, evolve, and ensure
`pandora` is behaving as we expect. If the quality of our test code decrease, the quality of our production code will
eventually decrease. So, make your tests readable, maintenable, clean, and without duplicated code.
These conventions help you to achieve this.

NB: this is a fairly incomplete and sometimes imprecise document, far from being good, it will
get better with time. Some stuff here might not be the best way to do something, we will be improving it as we write more
 and more unit tests following these conventions.
 The big majority of `pandora` tests do not follow these conventions, they will be slowly refactored.

## Conventions

If you have any suggestions for these conventions, please feel free to add.

* Use mocked objects to emulate complex dependencies which are hard to setup;
  * For simple dependencies (e.g. an `Interval` object), feel free to skip mocking, as it is probably easier and readable to setup such simple dependencies with real objects; 

* Use dependency injection when designing classes to make it easier to mock during tests;

* Do not refer directly to test code inside production code
  * For example, if you have to test a `private` method of your class, make it `virtual protected`, create a mock class
  that extends your class to be tested, and change the visibility of the tested method to `public`. Prefer this than making test cases `friend`s
  of the tested class. There is a specific place to create such a mock class (in a Fixture class, see below);

* For each tested method of a class, always declare a Fixture class even if your tests do not have share a common setup/teardown
  * This will encourage the non-duplication of setup and teardowns, and makes it easier to understand each test as they follow the same standard;
  * The Fixture class naming should obey the following convention: `<class_name>Test___<method_name>___Fixture`
  
* Your fixture class should mock all dependencies of the tested method
  * This includes methods of the same class of the tested method;
  * Your fixture class should have several mock classes defined inside it;
  
* Each test should follow this very simple algorithm, which is divided in four sections:
  1. Variable setups (aided by mocking and the Fixture class);
  2. Invocation of the tested method/function;
  3. Teardown (aided by the Fixture class)
  4. Assertion of the result;

* No section in each test should have more than a few lines. Sections 2 and 4 should preferably have a single line.
If your tests are long, there should be probably some code that could be duplicated between tests, or a hidden concept that could be better
described in a function or class to improve readability. Remember that your main task is to test the code, but also make it very easy for anyone reading the test to
understand what it does;

* A test tests only one condition and nothing else; 

* Your test naming should obey the following convention: `TEST_F(<FixtureClassName>, <short_description_of_the_test>___<expected result>)`


## Templates

### Fixture class template

```
class <Class_name>Test___<Method_name>___Fixture : public ::testing::Test {
protected:
    class <First_mocked_class> : public <Real_class_or_another_mocked_class> {
    public:
        <methods_to_mock>
    };


    void SetUp() override {
        // Code here will be called right before each test
    }

    void TearDown() override {
        // Code here will be called immediately after each test 
    }

    // Objects declared here can be used by all tests for this fixture
    // You should have at least one object for each mocked class
};
```

### Test using the fixture class
```
TEST_F(<FixtureClassName>, <short_description_of_the_test>___<expected result>) {
    <Section_1: additional_setups>

    <Section_2: invocation_of_the_method>

    <Section_3: additional_teardowns>

    <Section_4: assertion_of_the_result>
}
```


## Examples
### Fixture class
```
class LocalPRGTest___get_number_of_bases_in_local_path_before_a_given_position___Fixture : public ::testing::Test {
protected:
    class LocalPRGMockExposesTestedMethod : public LocalPRGMock {
    public:
        virtual uint32_t get_number_of_bases_in_local_path_before_a_given_position(const std::vector<LocalNodePtr> &local_path,
                                                                                   uint32_t position) const {
            return LocalPRGMock::get_number_of_bases_in_local_path_before_a_given_position(local_path, position);
        }
    };


    void SetUp() override {
        local_path_with_two_intervals.push_back(std::make_shared<LocalNode>("", Interval(10, 20), 0));
        local_path_with_two_intervals.push_back(std::make_shared<LocalNode>("", Interval(35, 45), 1));
    }

    void TearDown() override {
    }

    LocalPRGMockExposesTestedMethod local_prg_mock;
    std::vector<LocalNodePtr> empty_local_path;
    std::vector<LocalNodePtr> local_path_with_two_intervals;
};
```

### Test
```
TEST_F(LocalPRGTest___get_number_of_bases_in_local_path_before_a_given_position___Fixture, empty_local_path___returns_0) {
    uint32_t position = 30;

    uint32_t actual = local_prg_mock.get_number_of_bases_in_local_path_before_a_given_position(empty_local_path, position);

    uint32_t expected = 0;
    ASSERT_EQ(actual, expected);
}
```