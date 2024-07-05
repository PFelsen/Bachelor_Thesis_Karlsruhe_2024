#include "Point.hpp"
#include "TestEnvironment.hpp"
#include "utility/Config.hpp"

using namespace std;

template<typename T>
struct ConfigTestParam {
  string key;
  string valueString;
  T exptectedValue;
};

template<typename T>
class ConfigTest : public TestWithParam<ConfigTestParam<T>> {
protected:
  string key;
  string valueString;
  T expectedValue;

  ConfigTest(ConfigTestParam<T> param) :
      key(param.key), valueString(param.valueString), expectedValue(param.exptectedValue) {
    map<string, string> configMap{{key, valueString}};
    Config::Initialize(configMap);
  }

  void SetUp() override {}

  void TearDown() override { Config::Close(); }

  void checkExpectedValue() {
    T result = Config::GetWithoutCheck<T>(key);
    EXPECT_EQ(expectedValue, result) << "key: " + key;
  }

  void checkExpectedSizeForIterable() {
    T result = Config::GetWithoutCheck<T>(key);
    EXPECT_EQ(expectedValue.size(), result.size()) << "key: " + key;
  }
};

class ConfigTestBool : public ConfigTest<bool> {
public:
  ConfigTestBool() : ConfigTest(GetParam()){};
};

TEST_P(ConfigTestBool, TestBool) { checkExpectedValue(); }

INSTANTIATE_TEST_SUITE_P(ConfigTest, ConfigTestBool,
                         Values(ConfigTestParam<bool>{"myBool1", "1", true},
                                ConfigTestParam<bool>{"myBool2", "0", false},
                                ConfigTestParam<bool>{"myBool3", "true", true},
                                ConfigTestParam<bool>{"myBool4", "True", true},
                                ConfigTestParam<bool>{"myBool5", "false", false},
                                ConfigTestParam<bool>{"myBool6", "False", false}));

class ConfigTestVectorDouble : public ConfigTest<vector<double>> {
public:
  ConfigTestVectorDouble() : ConfigTest(GetParam()){};
};

TEST_P(ConfigTestVectorDouble, TestVector) {
  checkExpectedValue();
  checkExpectedSizeForIterable();
}

INSTANTIATE_TEST_SUITE_P(
    ConfigTest, ConfigTestVectorDouble,
    Values(ConfigTestParam<vector<double>>{"myVector1", "1.05,-1.00003,1334.5",
                                           vector<double>{1.05, -1.00003, 1334.5}},
           ConfigTestParam<vector<double>>{"myVector2", "1,2,3,4", vector<double>{1, 2, 3, 4}},
           ConfigTestParam<vector<double>>{"myVector3", "-1,-2222", vector<double>{-1, -2222}},
           ConfigTestParam<vector<double>>{"myVector4", "1.000043,456786.34534",
                                           vector<double>{1.000043, 456786.34534}},
           ConfigTestParam<vector<double>>{"myVector5", "[1.000043, 456786.34534]",
                                           vector<double>{1.000043, 456786.34534}}));

class ConfigTestVectorInt : public ConfigTest<vector<int>> {
public:
  ConfigTestVectorInt() : ConfigTest(GetParam()){};
};

TEST_P(ConfigTestVectorInt, TestVector) {
  checkExpectedValue();
  checkExpectedSizeForIterable();
}

INSTANTIATE_TEST_SUITE_P(
    ConfigTest, ConfigTestVectorInt,
    Values(ConfigTestParam<vector<int>>{"myVector1", "", vector<int>{}},
           ConfigTestParam<vector<int>>{"myVector2", "1", vector<int>{1}},
           ConfigTestParam<vector<int>>{"myVector3", "1,2,3,4", vector<int>{1, 2, 3, 4}},
           ConfigTestParam<vector<int>>{"myVector4", "1,2,3,4,", vector<int>{1, 2, 3, 4}},
           ConfigTestParam<vector<int>>{"myVector5", ",  ,1,   ,,2, ,,3,4, ,  ,,  ",
                                        vector<int>{1, 2, 3, 4}},
           ConfigTestParam<vector<int>>{"myVector6", "-1,-2,-3,-4", vector<int>{-1, -2, -3, -4}},
           ConfigTestParam<vector<int>>{"myVector7", "1,-1,2,-2,3,-3",
                                        vector<int>{1, -1, 2, -2, 3, -3}},
           ConfigTestParam<vector<int>>{"myVector8", "34622457,42437257,-24652577",
                                        vector<int>{34622457, 42437257, -24652577}},
           ConfigTestParam<vector<int>>{"myVector8", "[3,4,5]", vector<int>{3, 4, 5}}

           ));

class ConfigTestDouble : public ConfigTest<double> {
public:
  ConfigTestDouble() : ConfigTest(GetParam()){};
};

TEST_P(ConfigTestDouble, TestDouble) { checkExpectedValue(); }

INSTANTIATE_TEST_SUITE_P(ConfigTest, ConfigTestDouble,
                         Values(ConfigTestParam<double>{"myDouble1", "1.534536", 1.534536},
                                ConfigTestParam<double>{"myDouble2", "0", 0.0},
                                ConfigTestParam<double>{"myDouble3", "-1445.8989555",
                                                        -1445.8989555},
                                ConfigTestParam<double>{"myDouble4", "545", 545.0}));

class ConfigTestString : public ConfigTest<string> {
public:
  ConfigTestString() : ConfigTest(GetParam()){};
};

TEST_P(ConfigTestString, TestString) { checkExpectedValue(); }

INSTANTIATE_TEST_SUITE_P(ConfigTest, ConfigTestString,
                         Values(ConfigTestParam<string>{"myString1", "  White  Spaces  ",
                                                        "White  Spaces"},
                                ConfigTestParam<string>{"myString2", "     linear", "linear"},
                                ConfigTestParam<string>{"myString3", "1445     ", "1445"},
                                ConfigTestParam<string>{"myString4", "false", "false"}));

class ConfigTestPoint : public ConfigTest<Point> {
public:
  ConfigTestPoint() : ConfigTest(GetParam()){};

  void checkPoint() {
    Point result = Config::GetWithoutCheck<Point>(key);
    EXPECT_EQ(expectedValue, result) << "key: " + key;
  }
};

TEST_P(ConfigTestPoint, TestPoint) { checkExpectedValue(); }

INSTANTIATE_TEST_SUITE_P(ConfigTest, ConfigTestPoint,
                         Values(ConfigTestParam<Point>{"myPoint1", "1,2,3,4", Point(1, 2, 3, 4)},
                                ConfigTestParam<Point>{"myPoint2", "1,2,3", Point(1, 2, 3, 0)},
                                ConfigTestParam<Point>{"myPoint3", "1,2", Point(1, 2)},
                                ConfigTestParam<Point>{"myPoint4", "1", Point(1)}));

TEST(ConfigParseTest, ParseTest) {
  string conf = "a=b;c=d;";
  std::unordered_map<string, ConfigEntry> map = Config::toMap(conf);

  ASSERT_EQ(map.size(), 2);
  ASSERT_EQ(map.find("a")->second.value, "b");
  ASSERT_EQ(map.find("c")->second.value, "d");
}

TEST(ConfigParseTest, OrderParseTest) {
  string conf = "a=b;a=c;";
  std::unordered_map<string, ConfigEntry> map = Config::toMap(conf);

  ASSERT_EQ(map.size(), 1);
  ASSERT_EQ(map.find("a")->second.value, "b");
}

TEST(ConfigParseTest, ReferenceParseTest) {
  string conf = "a=b;c=&a;";
  std::unordered_map<string, ConfigEntry> map = Config::toMap(conf);
  Config::parseReferenceConfig(map);
  ASSERT_EQ(map.size(), 2);
  ASSERT_EQ(map.find("a")->second.value, "b");
  ASSERT_EQ(map.find("c")->second.value, "b");
}

int main(int argc, char **argv) {
  MppTest mppTest = MppTestBuilder(argc, argv).WithoutDefaultConfig();
  return mppTest.RUN_ALL_MPP_TESTS();
}
