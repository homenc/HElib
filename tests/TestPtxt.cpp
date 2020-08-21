/* Copyright (C) 2019-2020 IBM Corp.
 * This program is Licensed under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. See accompanying LICENSE file.
 */

#include <numeric>

#include <helib/Ptxt.h>
#include <helib/helib.h>
#include <helib/replicate.h>
#include <helib/NumbTh.h>

#include "test_common.h"
#include "gtest/gtest.h"

namespace {

struct BGVParameters
{
  BGVParameters(unsigned m, unsigned p, unsigned r) : m(m), p(p), r(r){};

  const unsigned m;
  const unsigned p;
  const unsigned r;

  friend std::ostream& operator<<(std::ostream& os, const BGVParameters& params)
  {
    return os << "{"
              << "m = " << params.m << ", "
              << "p = " << params.p << ", "
              << "r = " << params.r << "}";
  }
};

class TestPtxtCKKS : public ::testing::TestWithParam<unsigned>
{
protected:
  TestPtxtCKKS() :
      // Only relevant parameter is m for a CKKS plaintext
      m(GetParam()),
      context(m, -1, 50)
  {}

  const unsigned long m;

  helib::Context context;

  const double pre_encryption_epsilon = 1E-8;
  const double post_encryption_epsilon = 1E-3;
};

TEST_P(TestPtxtCKKS, canBeConstructedWithCKKSContext)
{
  helib::Ptxt<helib::CKKS> ptxt(context);
}

TEST_P(TestPtxtCKKS, canBeDefaultConstructed) { helib::Ptxt<helib::CKKS> ptxt; }

TEST_P(TestPtxtCKKS, canBeCopyConstructed)
{
  helib::Ptxt<helib::CKKS> ptxt(context);
  helib::Ptxt<helib::CKKS> ptxt2(ptxt);
}

TEST_P(TestPtxtCKKS, canBeAssignedFromOtherPtxt)
{
  helib::Ptxt<helib::CKKS> ptxt(context);
  helib::Ptxt<helib::CKKS> ptxt2 = ptxt;
}

TEST_P(TestPtxtCKKS, reportsWhetherItIsValid)
{
  helib::Ptxt<helib::CKKS> invalid_ptxt;
  helib::Ptxt<helib::CKKS> valid_ptxt(context);
  EXPECT_FALSE(invalid_ptxt.isValid());
  EXPECT_TRUE(valid_ptxt.isValid());
}

TEST_P(TestPtxtCKKS, hasSameNumberOfSlotsAsContext)
{
  helib::Ptxt<helib::CKKS> ptxt(context);
  EXPECT_EQ(ptxt.size(), context.ea->size());
}

TEST_P(TestPtxtCKKS, preservesDataPassedIntoConstructor)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (std::size_t i = 0; i < data.size(); ++i)
    data[i] = i / 10.0;
  helib::Ptxt<helib::CKKS> ptxt(context, data);

  EXPECT_EQ(ptxt.size(), data.size());
  for (std::size_t i = 0; i < data.size(); ++i)
    EXPECT_EQ(ptxt[i], data[i]);
}

TEST_P(TestPtxtCKKS, hasSameNumberOfSlotsAsContextWhenCreatedWithData)
{
  std::vector<std::complex<double>> data(context.ea->size() - 1);
  helib::Ptxt<helib::CKKS> ptxt(context, data);
  EXPECT_EQ(ptxt.size(), context.ea->size());
}

TEST_P(TestPtxtCKKS, replicateValueWhenPassingASingleSlotTypeNumber)
{
  std::complex<double> num = {1. / 10.0, 1. / 10.0};

  helib::Ptxt<helib::CKKS> ptxt(context, num);
  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], num);
  }
}

TEST_P(TestPtxtCKKS, replicateValueWhenPassingASingleNonSlotTypeNumber)
{
  double num = 1. / 10.0;

  helib::Ptxt<helib::CKKS> ptxt(context, num);
  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], num);
  }
}

TEST_P(TestPtxtCKKS, atMethodThrowsOrReturnsCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (std::size_t i = 0; i < data.size(); ++i)
    data[i] = i / 10.0;
  helib::Ptxt<helib::CKKS> ptxt(context, data);

  for (long i = -5; i < 0; ++i) {
    EXPECT_THROW(ptxt.at(i), helib::OutOfRangeError);
  }
  for (long i = 0; i < helib::lsize(data); ++i) {
    EXPECT_EQ(ptxt.at(i), data.at(i));
  }
  for (std::size_t i = data.size(); i < data.size() + 5; ++i) {
    EXPECT_THROW(ptxt.at(i), helib::OutOfRangeError);
  }
}

TEST_P(TestPtxtCKKS, padsWithZerosWhenPassingInSmallDataVector)
{
  std::vector<std::complex<double>> data(context.ea->size() - 1);
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = {(i - 1) / 10.0, (i - 1) / 10.0};
  }
  helib::Ptxt<helib::CKKS> ptxt(context, data);
  for (std::size_t i = 0; i < data.size(); ++i) {
    EXPECT_EQ(ptxt[i], data[i]);
  }
  for (std::size_t i = data.size(); i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], 0.0);
  }
}

TEST_P(TestPtxtCKKS, preservesDataPassedIntoConstructorAsDouble)
{
  std::vector<double> data(context.ea->size() - 1);
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = (i - 1) / 10.0;
  }
  helib::Ptxt<helib::CKKS> ptxt(context, data);
  for (std::size_t i = 0; i < data.size(); ++i) {
    EXPECT_EQ(ptxt[i], data[i]);
  }
  for (std::size_t i = data.size(); i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], 0.0);
  }
}

TEST_P(TestPtxtCKKS, writesDataCorrectlyToOstream)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {(i * i) / 10.0, (i * i * i) / 7.5};
  }
  helib::Ptxt<helib::CKKS> ptxt(context, data);
  std::stringstream ss;
  ss << "[";
  ss << std::setprecision(std::numeric_limits<double>::digits10);
  for (auto it = data.begin(); it != data.end(); it++) {
    ss << "[" << it->real() << ", " << it->imag() << "]";
    if (it != data.end() - 1) {
      ss << ", ";
    }
  }
  ss << "]";
  std::string expected = ss.str();
  std::ostringstream os;
  os << ptxt;

  EXPECT_EQ(os.str(), expected);
}

TEST_P(TestPtxtCKKS, readsDataCorrectlyFromIstream)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {(i * i) / 10.0, (i * i * i) / 7.5};
  }
  helib::Ptxt<helib::CKKS> ptxt(context);
  std::stringstream ss;
  ss << "[";
  ss << std::setprecision(std::numeric_limits<double>::digits10);
  for (auto it = data.begin(); it != data.end(); it++) {
    helib::serialize(ss, *it);
    if (it != data.end() - 1) {
      ss << ", ";
    }
  }
  ss << "]";
  std::istringstream is(ss.str());
  is >> ptxt;

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_NEAR(std::abs(ptxt[i] - data[i]), 0, pre_encryption_epsilon);
  }
}

TEST_P(TestPtxtCKKS, readsSquareBracketsDataCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {(i * i) / 10.0, (i * i * i) / 7.5};
  }
  helib::Ptxt<helib::CKKS> ptxt(context);
  std::stringstream ss;
  ss << "[";
  ss << std::setprecision(std::numeric_limits<double>::digits10);
  for (auto it = data.begin(); it != data.end(); it++) {
    ss << "[" << it->real() << ", " << it->imag() << "]";
    if (it != data.end() - 1) {
      ss << ", ";
    }
  }
  ss << "]";
  std::istringstream is(ss.str());
  is >> ptxt;

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_NEAR(std::abs(ptxt[i] - data[i]), 0, pre_encryption_epsilon);
  }
}

TEST_P(TestPtxtCKKS, serializeFunctionSerializesStdComplexCorrectly)
{
  // TODO: This test may be removed from the fixture and put as standalone
  std::stringstream ss;
  std::complex<double> num;

  num = 0.0;
  helib::serialize(ss, num);
  EXPECT_EQ(ss.str(), "[0, 0]");
  ss.str("");

  num = 10.3;
  helib::serialize(ss, num);
  EXPECT_EQ(ss.str(), "[10.3, 0]");
  ss.str("");

  num = {0, 10.3};
  helib::serialize(ss, num);
  EXPECT_EQ(ss.str(), "[0, 10.3]");
  ss.str("");

  num = {3.3, 16.6};
  helib::serialize(ss, num);
  EXPECT_EQ(ss.str(), "[3.3, 16.6]");
  ss.str("");
}

TEST_P(TestPtxtCKKS, deserializeFunctionDeserializesStdComplexCorrectly)
{
  // TODO: This test may be removed from the fixture and put as standalone
  std::complex<double> num, expected;
  std::stringstream ss;

  num = 0.0;
  ss << "[1,2,3]";
  expected = 0.0;
  EXPECT_THROW(helib::deserialize(ss, num), helib::IOError);
  ss.str("");

  num = 0.0;
  ss << "[]";
  expected = 0.0;
  helib::deserialize(ss, num);
  EXPECT_NEAR(std::abs(num - expected), 0.0, pre_encryption_epsilon);
  ss.str("");

  num = 0.0;
  ss << "[0.0]";
  expected = 0.0;
  helib::deserialize(ss, num);
  EXPECT_NEAR(std::abs(num - expected), 0.0, pre_encryption_epsilon);
  ss.str("");

  num = 0.0;
  ss << "[0.0,0.0]";
  expected = 0.0;
  helib::deserialize(ss, num);
  EXPECT_NEAR(std::abs(num - expected), 0.0, pre_encryption_epsilon);
  ss.str("");

  num = 0.0;
  ss << "[5.3,0]";
  expected = 5.3;
  helib::deserialize(ss, num);
  EXPECT_NEAR(std::abs(num - expected), 0.0, pre_encryption_epsilon);
  ss.str("");

  num = 0.0;
  ss << "[0,8.16]";
  expected = {0.0, 8.16};
  helib::deserialize(ss, num);
  EXPECT_NEAR(std::abs(num - expected), 0.0, pre_encryption_epsilon);
  ss.str("");

  num = 0.0;
  ss << "[3.4,9.99]";
  expected = {3.4, 9.99};
  helib::deserialize(ss, num);
  EXPECT_NEAR(std::abs(num - expected), 0.0, pre_encryption_epsilon);
  ss.str("");
}

TEST_P(TestPtxtCKKS, serializeFunctionSerializesCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {(i * i) / 10.0, (i * i * i) / 7.5};
  }
  helib::Ptxt<helib::CKKS> ptxt(context, data);
  std::stringstream ss;
  ss << "[";
  ss << std::setprecision(std::numeric_limits<double>::digits10);
  for (auto it = data.begin(); it != data.end(); it++) {
    ss << "[" << it->real() << ", " << it->imag() << "]";
    if (it != data.end() - 1) {
      ss << ", ";
    }
  }
  ss << "]";
  std::stringstream serialized_ptxt;
  helib::serialize(serialized_ptxt, ptxt);

  EXPECT_EQ(serialized_ptxt.str(), ss.str());
}

TEST_P(TestPtxtCKKS, deserializeWorksCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {(i * i) / 10.0, (i * i * i) / 7.5};
  }
  helib::Ptxt<helib::CKKS> ptxt(context);
  std::stringstream ss;
  ss << "[";
  ss << std::setprecision(std::numeric_limits<double>::digits10);
  for (auto it = data.begin(); it != data.end(); it++) {
    ss << "[" << it->real() << ", " << it->imag() << "]";
    if (it != data.end() - 1) {
      ss << ", ";
    }
  }
  ss << "]";
  std::istringstream is(ss.str());
  is >> ptxt;

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_NEAR(std::abs(ptxt[i] - data[i]), 0, pre_encryption_epsilon);
  }
}

TEST_P(TestPtxtCKKS, deserializeFunctionThrowsIfMoreElementsThanSlots)
{
  std::vector<std::complex<double>> data(context.ea->size() + 1);
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {(i * i) / 10.0, (i * i * i) / 7.5};
  }
  helib::Ptxt<helib::CKKS> ptxt(context);
  std::stringstream ss;
  ss << "[";
  ss << std::setprecision(std::numeric_limits<double>::digits10);
  for (auto it = data.begin(); it != data.end(); it++) {
    ss << "[" << it->real() << ", " << it->imag() << "]";
    if (it != data.end() - 1) {
      ss << ", ";
    }
  }
  ss << "]";
  std::istringstream is(ss.str());
  EXPECT_THROW(helib::deserialize(is, ptxt), helib::IOError);
}

TEST_P(TestPtxtCKKS, rightShiftOperatorThrowsIfMoreElementsThanSlots)
{
  std::vector<std::complex<double>> data(context.ea->size() + 1);
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {(i * i) / 10.0, (i * i * i) / 7.5};
  }
  helib::Ptxt<helib::CKKS> ptxt(context);
  std::stringstream ss;
  ss << "[";
  ss << std::setprecision(std::numeric_limits<double>::digits10);
  for (auto it = data.begin(); it != data.end(); it++) {
    ss << "[" << it->real() << ", " << it->imag() << "]";
    if (it != data.end() - 1) {
      ss << ", ";
    }
  }
  ss << "]";
  std::istringstream is(ss.str());
  EXPECT_THROW(is >> ptxt, helib::IOError);
}

TEST_P(TestPtxtCKKS, deserializeIsInverseOfSerialize)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {(i * i) / 10.0, (i * i * i) / 7.5};
  }
  helib::Ptxt<helib::CKKS> ptxt(context, data);

  std::stringstream str;
  str << ptxt;

  helib::Ptxt<helib::CKKS> deserialized(context);
  str >> deserialized;

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_NEAR(std::abs(ptxt[i] - deserialized[i]), 0, pre_encryption_epsilon);
  }
}

TEST_P(TestPtxtCKKS, readsManyPtxtsFromStream)
{
  std::vector<std::complex<double>> data1(context.ea->size());
  std::vector<std::complex<double>> data2(context.ea->size());
  std::vector<std::complex<double>> data3(context.ea->size());
  for (std::size_t i = 0; i < data1.size(); ++i) {
    data1[i] = {(i * i) / 10.0, (i * i * i) / 7.5};
    data2[i] = data1[i] * 2.0;
    data3[i] = data1[i] * 3.5;
  }
  helib::Ptxt<helib::CKKS> ptxt1(context, data1);
  helib::Ptxt<helib::CKKS> ptxt2(context, data2);
  helib::Ptxt<helib::CKKS> ptxt3(context, data3);

  std::stringstream ss;
  ss << ptxt1 << std::endl;
  ss << ptxt2 << std::endl;
  ss << ptxt3 << std::endl;

  helib::Ptxt<helib::CKKS> deserialized1(context);
  helib::Ptxt<helib::CKKS> deserialized2(context);
  helib::Ptxt<helib::CKKS> deserialized3(context);
  ss >> deserialized1;
  ss >> deserialized2;
  ss >> deserialized3;

  for (std::size_t i = 0; i < ptxt1.size(); ++i) {
    EXPECT_NEAR(std::abs(ptxt1[i] - deserialized1[i]),
                0,
                pre_encryption_epsilon);
    EXPECT_NEAR(std::abs(ptxt2[i] - deserialized2[i]),
                0,
                pre_encryption_epsilon);
    EXPECT_NEAR(std::abs(ptxt3[i] - deserialized3[i]),
                0,
                pre_encryption_epsilon);
  }
}

TEST_P(TestPtxtCKKS, getSlotReprReturnsData)
{
  std::vector<std::complex<double>> data(context.ea->size() - 1);
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = {(i - 1) / 10.0, (i - 1) / 10.0};
  }
  helib::Ptxt<helib::CKKS> ptxt(context, data);
  std::vector<std::complex<double>> expected_repr(context.ea->size());
  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    expected_repr[i] = i < data.size() ? data[i] : 0;
  }
  EXPECT_EQ(ptxt.getSlotRepr(), expected_repr);
}

TEST_P(TestPtxtCKKS, runningSumsWorksCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {i / 1.0, (i * i) / 1.0};
  }

  helib::Ptxt<helib::CKKS> ptxt(context, data);
  ptxt.runningSums();

  std::vector<std::complex<double>> expected_result(context.ea->size());
  for (std::size_t i = 0; i < data.size(); ++i)
    expected_result[i] = {(i * (i + 1)) / 2.0,
                          (i * (i + 1) * (2 * i + 1)) / 6.0};

  EXPECT_EQ(ptxt.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, totalSumsWorksCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {i / 1.0, (i * i) / 1.0};
  }

  helib::Ptxt<helib::CKKS> ptxt(context, data);
  ptxt.totalSums();

  std::vector<std::complex<double>> expected_result(context.ea->size());
  for (std::size_t i = 0; i < data.size(); ++i)
    expected_result[i] = {
        ((data.size() - 1) * data.size()) / 2.0,
        ((data.size() - 1) * data.size() * (2 * (data.size() - 1) + 1)) / 6.0};

  EXPECT_EQ(ptxt.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, incrementalProductWorksCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {(i + 1) / 5.0, (i * i + 1) / 10.0};
  }

  helib::Ptxt<helib::CKKS> ptxt(context, data);
  ptxt.incrementalProduct();

  std::vector<std::complex<double>> expected_result(data);
  for (std::size_t i = 1; i < data.size(); ++i)
    expected_result[i] *= expected_result[i - 1];

  EXPECT_EQ(ptxt.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, totalProductWorksCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {(i + 1) / 10.0, (i * i + 1) / 10.0};
  }

  helib::Ptxt<helib::CKKS> ptxt(context, data);
  ptxt.totalProduct();

  std::complex<double> product = {1.0, 0.0};
  for (std::size_t i = 0; i < data.size(); ++i)
    product *= data[i];
  std::vector<std::complex<double>> expected_result(context.ea->size(),
                                                    product);

  EXPECT_EQ(ptxt.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, innerProductWorksCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {i / 1.0, (i * i) / 1.0};
  }

  helib::Ptxt<helib::CKKS> ptxt(context, data);
  std::vector<helib::Ptxt<helib::CKKS>> first_ptxt_vector(2, ptxt);
  ptxt += ptxt;
  std::vector<helib::Ptxt<helib::CKKS>> second_ptxt_vector(3, ptxt);

  helib::Ptxt<helib::CKKS> result(context);
  innerProduct(result, first_ptxt_vector, second_ptxt_vector);

  std::vector<std::complex<double>> expected_result(data.size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    expected_result[i] = (data[i] * (data[i] + data[i]));
    expected_result[i] +=
        expected_result[i]; // expected_result = 2*expected_result
  }

  EXPECT_EQ(result.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, mapTo01MapsSlotsCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {i / 1.0, (i * i) / 1.0};
  }

  helib::Ptxt<helib::CKKS> ptxt(context, data);
  helib::Ptxt<helib::CKKS> ptxt2(context, data);
  // Should exist as a free function and a member function
  ptxt.mapTo01();
  mapTo01(*(context.ea), ptxt2);

  std::vector<std::complex<double>> expected_result(context.ea->size());
  for (std::size_t i = 1; i < data.size(); ++i) {
    expected_result[i] = {1, 0};
  }

  EXPECT_EQ(ptxt.getSlotRepr(), expected_result);
  EXPECT_EQ(ptxt2.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, timesEqualsOtherPlaintextWorks)
{
  std::vector<std::complex<double>> product_data(context.ea->size(),
                                                 {-3.14, -1.0});
  std::vector<std::complex<double>> multiplier_data(context.ea->size());
  for (long i = 0; i < helib::lsize(multiplier_data); ++i) {
    multiplier_data[i] = {(i - 1) / 10.0, (i + 1) / 10.0};
  }

  std::vector<std::complex<double>> expected_result(product_data);
  for (std::size_t i = 0; i < product_data.size(); ++i) {
    expected_result[i] = expected_result[i] * multiplier_data[i];
  }

  helib::Ptxt<helib::CKKS> product(context, product_data);
  helib::Ptxt<helib::CKKS> multiplier(context, multiplier_data);

  product *= multiplier;

  EXPECT_EQ(product.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, minusEqualsOtherPlaintextWorks)
{
  std::vector<std::complex<double>> difference_data(context.ea->size(),
                                                    {2.718, -1.0});
  std::vector<std::complex<double>> subtrahend_data(context.ea->size());
  for (long i = 0; i < helib::lsize(subtrahend_data); ++i) {
    subtrahend_data[i] = {(i - 1) / 10.0, (i + 1) / 10.0};
  }

  std::vector<std::complex<double>> expected_result(difference_data);
  for (std::size_t i = 0; i < subtrahend_data.size(); ++i) {
    expected_result[i] = expected_result[i] - subtrahend_data[i];
  }

  helib::Ptxt<helib::CKKS> difference(context, difference_data);
  helib::Ptxt<helib::CKKS> subtrahend(context, subtrahend_data);

  difference -= subtrahend;

  EXPECT_EQ(difference.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, minusEqualsComplexScalarWorks)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = {(i * i - 1) / 10.0, (i * i + 1) / 10.0};
  }

  const std::complex<double> scalar = {2.5, -0.5};

  std::vector<std::complex<double>> expected_result(data);
  for (auto& num : expected_result)
    num = num - scalar;

  helib::Ptxt<helib::CKKS> ptxt(context, data);
  ptxt -= scalar;

  EXPECT_EQ(ptxt.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, minusEqualsNonComplexScalarWorks)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = {(i * i - 1) / 2.0, (i * i + 1) / 5.0};
  }

  const double scalar = 15.3;
  const int int_scalar = 2;

  std::vector<std::complex<double>> expected_result(data);
  for (auto& num : expected_result) {
    num = num - scalar - static_cast<double>(int_scalar);
  }

  helib::Ptxt<helib::CKKS> ptxt(context, data);
  ptxt -= scalar;
  ptxt -= int_scalar;

  EXPECT_EQ(ptxt.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, plusEqualsOtherPlaintextWorks)
{
  std::vector<std::complex<double>> augend_data(context.ea->size());
  std::vector<std::complex<double>> addend_data(context.ea->size());
  for (long i = 0; i < helib::lsize(addend_data); ++i) {
    augend_data[i] = {i / 10.0, i * i / 10.0};
    addend_data[i] = {i / 20.0, i * i / 20.0};
  }
  std::vector<std::complex<double>> expected_result(context.ea->size());
  for (std::size_t i = 0; i < expected_result.size(); ++i)
    expected_result[i] = augend_data[i] + addend_data[i];

  helib::Ptxt<helib::CKKS> sum(context, augend_data);
  helib::Ptxt<helib::CKKS> addend(context, addend_data);
  sum += addend;

  EXPECT_EQ(sum.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, plusEqualsComplexScalarWorks)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = {-i / 10.0, (3 - i) / 4.0};
  }

  const double scalar = 3.14;

  std::vector<std::complex<double>> expected_result(data);
  for (auto& num : expected_result)
    num += scalar;

  helib::Ptxt<helib::CKKS> ptxt(context, data);
  ptxt += scalar;

  EXPECT_EQ(ptxt.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, plusEqualsNonComplexScalarWorks)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = {i + i / 5.0, i - i / 4.0};
  }

  const double scalar = 3.28;
  const int int_scalar = 13;

  std::vector<std::complex<double>> expected_result(data);
  for (auto& num : expected_result)
    num += scalar + static_cast<double>(int_scalar);

  helib::Ptxt<helib::CKKS> ptxt(context, data);
  ptxt += scalar;
  ptxt += int_scalar;

  const auto& slots = ptxt.getSlotRepr();
  for (std::size_t i = 0; i < data.size(); ++i) {

    EXPECT_DOUBLE_EQ(slots[i].real(), expected_result[i].real());
    EXPECT_DOUBLE_EQ(slots[i].imag(), expected_result[i].imag());
  }
}

TEST_P(TestPtxtCKKS, timesEqualsScalarWorks)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = {i * i * i / 100.0, -i / 3.0};
  }

  const double scalar = 10.28;
  const int int_scalar = -2;

  std::vector<std::complex<double>> expected_result(data);
  for (auto& num : expected_result)
    num *= scalar * static_cast<double>(int_scalar);

  helib::Ptxt<helib::CKKS> ptxt(context, data);
  ptxt *= scalar;
  ptxt *= int_scalar;

  const auto& slots = ptxt.getSlotRepr();
  for (std::size_t i = 0; i < data.size(); ++i) {

    EXPECT_DOUBLE_EQ(slots[i].real(), expected_result[i].real());
    EXPECT_DOUBLE_EQ(slots[i].imag(), expected_result[i].imag());
  }
}

TEST_P(TestPtxtCKKS, equalityWithOtherPlaintextWorks)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = {i * 2.5, (i - 2) * 2.5};
  }
  helib::Ptxt<helib::CKKS> ptxt1(context, data);
  helib::Ptxt<helib::CKKS> ptxt2(context, data);
  EXPECT_TRUE(ptxt1 == ptxt2);
  EXPECT_FALSE(ptxt1 == helib::Ptxt<helib::CKKS>());
}

TEST_P(TestPtxtCKKS, notEqualsOperatorWithOtherPlaintextWorks)
{
  std::vector<std::complex<double>> data1(context.ea->size());
  std::vector<std::complex<double>> data2(context.ea->size());
  for (long i = 0; i < helib::lsize(data1); ++i) {
    data1[i] = {(i + 1) * 2.5,
                -i * 2.5}; // i+1 makes the first element differ from (0,0)
    data2[i] = {i * 2.5, i * 6.5};
  }
  helib::Ptxt<helib::CKKS> ptxt1(context, data1);
  helib::Ptxt<helib::CKKS> ptxt2(context, data2);
  EXPECT_TRUE(ptxt1 != ptxt2);
  EXPECT_FALSE(ptxt1 != ptxt1);
}

TEST_P(TestPtxtCKKS, negateNegatesCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  const double pi = std::acos(-1);
  for (long j = 0; j < helib::lsize(data); ++j) {
    // Spiral with j -> j e^{2*i*pi*j/data.size()}
    data[j] = std::complex<double>{static_cast<double>(j), 0} *
              std::exp(std::complex<double>{0, 2.0 * pi * j / data.size()});
  }

  std::vector<std::complex<double>> expected_result(data);
  for (auto& num : expected_result)
    num *= std::complex<double>{-1.0, 0};

  helib::Ptxt<helib::CKKS> ptxt(context, data);
  ptxt.negate();

  EXPECT_EQ(ptxt.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, addConstantWorksCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (long i = 0; i < helib::lsize(data); ++i)
    data[i] = {i * 4.5, i / 2.0};

  std::vector<std::complex<double>> expected_result(data);
  for (auto& num : expected_result)
    (num += 5) += std::complex<double>{0, 0.5};

  helib::Ptxt<helib::CKKS> ptxt(context, data);
  ptxt.addConstantCKKS(5).addConstantCKKS(std::complex<double>{0, 0.5});

  EXPECT_EQ(ptxt.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, multiplyByMultipliesCorrectly)
{
  std::vector<std::complex<double>> product_data(context.ea->size());
  std::vector<std::complex<double>> multiplier_data(context.ea->size());
  for (long i = 0; i < helib::lsize(multiplier_data); ++i) {
    product_data[i] = {(2 - i) / 10.0, (1 - i) / 10.0};
    multiplier_data[i] = {std::exp(i / 100.), std::cos(i) * 12};
  }

  std::vector<std::complex<double>> expected_result(product_data);
  for (std::size_t i = 0; i < product_data.size(); ++i) {
    expected_result[i] = expected_result[i] * multiplier_data[i];
  }

  helib::Ptxt<helib::CKKS> product(context, product_data);
  helib::Ptxt<helib::CKKS> multiplier(context, multiplier_data);

  product.multiplyBy(multiplier);

  EXPECT_EQ(product.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, multiplyBy2MultipliesCorrectly)
{
  std::vector<std::complex<double>> product_data(context.ea->size());
  std::vector<std::complex<double>> multiplier_data1(context.ea->size());
  std::vector<std::complex<double>> multiplier_data2(context.ea->size());
  for (long i = 0; i < helib::lsize(multiplier_data1); ++i) {
    product_data[i] = static_cast<double>(i) *
                      std::exp(std::complex<double>{0, static_cast<double>(i)});
    multiplier_data2[i] =
        static_cast<double>(i) *
        std::exp(std::complex<double>{0, static_cast<double>(-i)});
    ;
    multiplier_data2[i] =
        5.0 * std::exp(std::complex<double>{0, static_cast<double>(i)});
    ;
    ;
  }

  std::vector<std::complex<double>> expected_result(product_data);
  for (std::size_t i = 0; i < product_data.size(); ++i) {
    expected_result[i] =
        expected_result[i] * multiplier_data1[i] * multiplier_data2[i];
  }

  helib::Ptxt<helib::CKKS> product(context, product_data);
  helib::Ptxt<helib::CKKS> multiplier1(context, multiplier_data1);
  helib::Ptxt<helib::CKKS> multiplier2(context, multiplier_data2);

  product.multiplyBy2(multiplier1, multiplier2);

  EXPECT_EQ(product.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, squareSquaresCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (long i = 0; i < helib::lsize(data); ++i) {
    // Lemniscate of Bernoulli
    double theta = 2. * std::acos(-1) * i / data.size();
    data[i] = std::cos(2. * theta) * std::exp(std::complex<double>{0, theta});
  }
  std::vector<std::complex<double>> expected_result(data);
  for (auto& num : expected_result)
    num *= num;
  helib::Ptxt<helib::CKKS> ptxt(context, data);
  ptxt.square();
  EXPECT_EQ(ptxt.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, cubeCubesCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (long i = 0; i < helib::lsize(data); ++i) {
    // Catenary
    data[i] = {static_cast<double>(1. * i - data.size() / 2) / data.size(),
               std::cosh((i - data.size() / 2.) / data.size())};
  }
  std::vector<std::complex<double>> expected_result(data);
  for (auto& num : expected_result)
    num = num * num * num;
  helib::Ptxt<helib::CKKS> ptxt(context, data);
  ptxt.cube();
  EXPECT_EQ(ptxt.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, powerCorrectlyRaisesToPowers)
{
  std::vector<std::complex<double>> data(context.ea->size());
  const double pi = std::acos(-1);
  // Spiral inside the unit disk
  for (long j = 0; j < helib::lsize(data); ++j) {
    data[j] = std::complex<double>{j / (double)data.size()} *
              std::exp(std::complex<double>{0, 2.0 * pi * j / data.size()});
  }
  std::vector<long> exponents{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1500};

  const auto naive_power = [](std::complex<double> base,
                              unsigned long exponent) {
    if (exponent == 0)
      return std::complex<double>{1.0};
    auto result = base;
    while (--exponent)
      result *= base;
    return result;
  };

  for (const auto& exponent : exponents) {
    std::vector<std::complex<double>> expected_result(data);
    for (auto& num : expected_result)
      num = naive_power(num, exponent);
    helib::Ptxt<helib::CKKS> ptxt(context, data);
    ptxt.power(exponent);
    for (std::size_t i = 0; i < ptxt.size(); ++i) {
      EXPECT_NEAR(std::norm(ptxt[i] - expected_result[i]),
                  0,
                  pre_encryption_epsilon);
    }
  }

  // Make sure raising to 0 throws
  helib::Ptxt<helib::CKKS> ptxt(context, data);
  EXPECT_THROW(ptxt.power(0l), helib::InvalidArgument);
}

TEST_P(TestPtxtCKKS, shiftShiftsRightCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  std::vector<std::complex<double>> right_shifted_data(context.ea->size());
  const auto non_neg_mod = [](int x, int mod) {
    return ((x % mod) + mod) % mod;
  };
  for (int i = 0; i < helib::lsize(data); ++i) {
    if (i > 3) {
      right_shifted_data[i] = {
          static_cast<double>(non_neg_mod(i - 3, data.size())),
          static_cast<double>(non_neg_mod(i - 3, data.size()))};
    }
    data[i] = {static_cast<double>(i), static_cast<double>(i)};
  }
  helib::Ptxt<helib::CKKS> ptxt(context, data);

  ptxt.shift(3);
  EXPECT_EQ(ptxt.getSlotRepr(), right_shifted_data);
}

TEST_P(TestPtxtCKKS, shiftShiftsLeftCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  std::vector<std::complex<double>> left_shifted_data(context.ea->size());
  const auto non_neg_mod = [](int x, int mod) {
    return ((x % mod) + mod) % mod;
  };
  for (int i = 0; i < helib::lsize(data); ++i) {
    if (i < long(data.size()) - 3 && data.size() > 3) {
      left_shifted_data[i] = {
          static_cast<double>(non_neg_mod(i + 3, data.size())),
          static_cast<double>(non_neg_mod(i + 3, data.size()))};
    }
    data[i] = {static_cast<double>(i), static_cast<double>(i)};
  }
  helib::Ptxt<helib::CKKS> ptxt(context, data);

  ptxt.shift(-3);
  EXPECT_EQ(ptxt.getSlotRepr(), left_shifted_data);
}

TEST_P(TestPtxtCKKS, shift1DShiftsRightCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  std::vector<std::complex<double>> right_shifted_data(context.ea->size());
  const auto non_neg_mod = [](int x, int mod) {
    return ((x % mod) + mod) % mod;
  };
  for (int i = 0; i < helib::lsize(data); ++i) {
    if (i > 3) {
      right_shifted_data[i] = {
          static_cast<double>(non_neg_mod(i - 3, data.size())),
          static_cast<double>(non_neg_mod(i - 3, data.size()))};
    }
    data[i] = {static_cast<double>(i), static_cast<double>(i)};
  }
  helib::Ptxt<helib::CKKS> ptxt(context, data);

  ptxt.shift1D(0, 3);
  EXPECT_EQ(ptxt.getSlotRepr(), right_shifted_data);
}

TEST_P(TestPtxtCKKS, shift1DShiftsLeftCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  std::vector<std::complex<double>> left_shifted_data(context.ea->size());
  const auto non_neg_mod = [](int x, int mod) {
    return ((x % mod) + mod) % mod;
  };
  for (int i = 0; i < helib::lsize(data); ++i) {
    if (i < long(data.size()) - 3 && data.size() > 3) {
      left_shifted_data[i] = {
          static_cast<double>(non_neg_mod(i + 3, data.size())),
          static_cast<double>(non_neg_mod(i + 3, data.size()))};
    }
    data[i] = {static_cast<double>(i), static_cast<double>(i)};
  }
  helib::Ptxt<helib::CKKS> ptxt(context, data);

  ptxt.shift1D(0, -3);
  EXPECT_EQ(ptxt.getSlotRepr(), left_shifted_data);
}

// These tests are disabled since the methods are private.
// These can be useful if tweaking the logic of this area.
// TEST(TestPtxtBGV, coord_To_Index_Works)
// {
//   helib::Context context(32109, 4999, 1);
//   helib::Ptxt<helib::BGV> ptxt(context);
//   std::vector<long> indices;
//   for(long i=0; i<6; ++i)
//     for(long j=0; j<2; ++j)
//       for(long k=0; k<2; ++k)
//         indices.push_back(ptxt.coordToIndex({i,j,k}));
//   std::vector<long> expected_indices(context.ea->size());
//   std::iota(expected_indices.begin(), expected_indices.end(), 0);
//   EXPECT_EQ(expected_indices, indices);
// }
//
// TEST(TestPtxtBGV, index_To_Coord_Works)
// {
//   helib::Context context(32109, 4999, 1);
//   helib::Ptxt<helib::BGV> ptxt(context);
//   std::vector<std::vector<long>> coords;
//   for(std::size_t i=0; i<ptxt.size(); ++i)
//     coords.push_back(ptxt.indexToCoord(i));
//   std::vector<std::vector<long>> expected_coords;
//   for(long i=0; i<6; ++i)
//     for(long j=0; j<2; ++j)
//       for(long k=0; k<2; ++k)
//         expected_coords.push_back({i,j,k});
//   EXPECT_EQ(expected_coords, coords);
// }

TEST_P(TestPtxtCKKS, rotate1DRotatesCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  std::vector<std::complex<double>> left_rotated_data(context.ea->size());
  const auto non_neg_mod = [](int x, int mod) {
    return ((x % mod) + mod) % mod;
  };
  for (int i = 0; i < helib::lsize(data); ++i) {
    data[i] = {static_cast<double>(non_neg_mod(i - 3, data.size())),
               static_cast<double>(non_neg_mod(i - 3, data.size()))};
    left_rotated_data[i] = {static_cast<double>(i), static_cast<double>(i)};
  }
  helib::Ptxt<helib::CKKS> ptxt(context, data);

  ptxt.rotate1D(0, -3);
  EXPECT_EQ(ptxt.getSlotRepr(), left_rotated_data);
  ptxt.rotate1D(0, 3);
  // Rotating back and forth gives the original data back
  EXPECT_EQ(ptxt.getSlotRepr(), data);
}

TEST_P(TestPtxtCKKS, rotateRotatesCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  std::vector<std::complex<double>> left_rotated_data(context.ea->size());
  const auto non_neg_mod = [](int x, int mod) {
    return ((x % mod) + mod) % mod;
  };
  for (int i = 0; i < helib::lsize(data); ++i) {
    data[i] = {static_cast<double>(non_neg_mod(i - 3, data.size())),
               static_cast<double>(non_neg_mod(i - 3, data.size()))};
    left_rotated_data[i] = {static_cast<double>(i), static_cast<double>(i)};
  }
  helib::Ptxt<helib::CKKS> ptxt(context, data);

  ptxt.rotate(-3);
  EXPECT_EQ(ptxt.getSlotRepr(), left_rotated_data);
  ptxt.rotate(3);
  // Rotating back and forth gives the original data back
  EXPECT_EQ(ptxt.getSlotRepr(), data);
}

TEST_P(TestPtxtCKKS, automorphWorksCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  std::vector<std::complex<double>> left_rotated_data(context.ea->size());
  const auto non_neg_mod = [](int x, int mod) {
    return ((x % mod) + mod) % mod;
  };
  for (int i = 0; i < helib::lsize(data); ++i) {
    data[i] = {static_cast<double>(non_neg_mod(i - 3, data.size())),
               static_cast<double>(non_neg_mod(i - 3, data.size()))};
  }
  helib::Ptxt<helib::CKKS> ptxt(context, data);
  helib::Ptxt<helib::CKKS> expected_result(context, data);

  long k = context.zMStar.ith_rep(1) ? context.zMStar.ith_rep(1) : 1;
  ptxt.automorph(k);
  expected_result.rotate(1);
  EXPECT_EQ(ptxt, expected_result);

  ptxt.automorph(context.zMStar.ith_rep(context.ea->size() - 1));
  expected_result.rotate(-1);
  EXPECT_EQ(ptxt, expected_result);
}

TEST_P(TestPtxtCKKS, replicateReplicatesCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = {i / 10.0, -i / 20.0};
  }
  helib::Ptxt<helib::CKKS> ptxt(context, data);
  helib::replicate(*context.ea, ptxt, data.size() - 1);
  std::vector<std::complex<double>> replicated_data(context.ea->size(),
                                                    data[data.size() - 1]);
  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], replicated_data[i]);
  }
}

TEST_P(TestPtxtCKKS, replicateAllWorksCorrectly)
{
  std::vector<std::complex<double>> data(context.ea->size());
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = {i / 10.0, -i / 20.0};
  }
  helib::Ptxt<helib::CKKS> ptxt(context, data);
  std::vector<helib::Ptxt<helib::CKKS>> replicated_ptxts = ptxt.replicateAll();
  for (long i = 0; i < helib::lsize(data); ++i) {
    for (const auto& slot : replicated_ptxts[i].getSlotRepr()) {
      EXPECT_EQ(data[i], slot);
    }
  }
}

TEST_P(TestPtxtCKKS, randomSetsDataRandomly)
{
  helib::Ptxt<helib::CKKS> ptxt(context);
  ptxt.random();
  std::vector<helib::Ptxt<helib::CKKS>> ptxts(5, ptxt);
  for (auto& p : ptxts)
    p.random();

  bool all_equal = true;
  for (std::size_t i = 0; i < ptxts.size() - 1; ++i)
    if (ptxts[i] != ptxts[i + 1]) {
      all_equal = false;
      break;
    }
  EXPECT_FALSE(all_equal) << "5 random ptxts are all equal - likely that"
                             " random() is not actually randomising!";
}

TEST_P(TestPtxtCKKS, complexConjCorrectlyConjugates)
{
  std::vector<std::complex<double>> data(context.ea->size());
  const std::complex<double> z{1, -1};
  for (long j = 0; j < helib::lsize(data); ++j) {
    // Line segment starting at 1 - i with gradient 2
    data[j] = z + std::complex<double>{j / 4.0, j / 2.0};
  }

  std::vector<std::complex<double>> expected_result(data);

  helib::Ptxt<helib::CKKS> ptxt(context, data);
  ptxt.complexConj();

  for (auto& num : expected_result)
    num = std::conj(num);

  EXPECT_EQ(ptxt.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, extractRealPartIsCorrect)
{
  std::vector<std::complex<double>> data(context.ea->size());
  const std::complex<double> z{1, -1};
  for (long j = 0; j < helib::lsize(data); ++j) {
    // Line segment starting at 1 - i with gradient 2
    data[j] = z + std::complex<double>{j / 4.0, j / 2.0};
  }

  std::vector<std::complex<double>> expected_result(data);
  for (auto& num : expected_result)
    num = std::real(num);

  helib::Ptxt<helib::CKKS> ptxt(context, data);
  context.ea->getCx().extractRealPart(ptxt);

  EXPECT_EQ(ptxt.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, extractImPartIsCorrect)
{
  std::vector<std::complex<double>> data(context.ea->size());
  const std::complex<double> z{1, -1};
  for (long j = 0; j < helib::lsize(data); ++j) {
    // Line segment starting at 1 - i with gradient 2
    data[j] = z + std::complex<double>{j / 4.0, j / 2.0};
  }

  std::vector<std::complex<double>> expected_result(data);
  for (auto& num : expected_result)
    num = std::imag(num);

  helib::Ptxt<helib::CKKS> ptxt(context, data);
  context.ea->getCx().extractImPart(ptxt);

  EXPECT_EQ(ptxt.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, realExtractsRealPart)
{
  std::vector<std::complex<double>> data(context.ea->size());
  const std::complex<double> z{1, -1};
  for (long j = 0; j < helib::lsize(data); ++j) {
    // Line segment starting at 1 - i with gradient 2
    data[j] = z + std::complex<double>{j / 4.0, j / 2.0};
  }

  std::vector<std::complex<double>> expected_result(data);
  for (auto& num : expected_result)
    num = std::real(num);

  helib::Ptxt<helib::CKKS> ptxt(context, data);
  ptxt = ptxt.real();

  EXPECT_EQ(ptxt.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, imagExtractsImaginaryPart)
{
  std::vector<std::complex<double>> data(context.ea->size());
  const std::complex<double> z{1, -1};
  for (long j = 0; j < helib::lsize(data); ++j) {
    // Line segment starting at 1 - i with gradient 2
    data[j] = z + std::complex<double>{j / 4.0, j / 2.0};
  }

  std::vector<std::complex<double>> expected_result(data);
  for (auto& num : expected_result)
    num = std::imag(num);

  helib::Ptxt<helib::CKKS> ptxt(context, data);
  ptxt = ptxt.imag();

  EXPECT_EQ(ptxt.getSlotRepr(), expected_result);
}

TEST_P(TestPtxtCKKS, canEncryptAndDecryptComplexPtxtsWithKeys)
{
  helib::buildModChain(context, 100, 2);
  helib::SecKey secret_key(context);
  secret_key.GenSecKey();
  const helib::PubKey& public_key(secret_key);

  std::vector<std::complex<double>> data(context.ea->size());
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = {(i - 3) / 5.0, i + 10.0};
  }
  helib::Ptxt<helib::CKKS> pre_encryption(context, data);
  helib::Ctxt ctxt(public_key);

  public_key.Encrypt(ctxt, pre_encryption);

  helib::Ptxt<helib::CKKS> post_decryption(context);
  secret_key.Decrypt(post_decryption, ctxt);
  EXPECT_EQ(pre_encryption.size(), post_decryption.size());
  for (std::size_t i = 0; i < pre_encryption.size(); ++i) {
    EXPECT_NEAR(std::norm(pre_encryption[i] - post_decryption[i]),
                0,
                post_encryption_epsilon);
  }
}

TEST_P(TestPtxtCKKS, canEncryptAndDecryptRealPtxtsWithKeys)
{
  helib::buildModChain(context, 100, 2);
  helib::SecKey secret_key(context);
  secret_key.GenSecKey();
  const helib::PubKey& public_key(secret_key);

  std::vector<double> data(context.ea->size());
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = (i - 3) / 10.0;
  }
  helib::Ptxt<helib::CKKS> pre_encryption(context, data);
  helib::Ctxt ctxt(public_key);

  public_key.Encrypt(ctxt, pre_encryption);

  helib::Ptxt<helib::CKKS> post_decryption(context);
  secret_key.Decrypt(post_decryption, ctxt);
  EXPECT_EQ(pre_encryption.size(), post_decryption.size());
  for (std::size_t i = 0; i < pre_encryption.size(); ++i) {
    EXPECT_NEAR(pre_encryption[i].real(),
                post_decryption[i].real(),
                post_encryption_epsilon);
    EXPECT_NEAR(pre_encryption[i].imag(),
                post_decryption[i].imag(),
                post_encryption_epsilon);
  }
}

TEST_P(TestPtxtCKKS, canEncryptAndDecryptComplexPtxtsWithEa)
{
  helib::buildModChain(context, 100, 2);
  helib::SecKey secret_key(context);
  secret_key.GenSecKey();
  const helib::PubKey& public_key(secret_key);

  std::vector<std::complex<double>> data(context.ea->size());
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = {(i - 3) / 10.0, i + 5.0};
  }
  helib::Ptxt<helib::CKKS> pre_encryption(context, data);
  helib::Ctxt ctxt(public_key);

  public_key.Encrypt(ctxt, pre_encryption);

  helib::Ptxt<helib::CKKS> post_decryption(context);
  secret_key.Decrypt(post_decryption, ctxt);
  EXPECT_EQ(pre_encryption.size(), post_decryption.size());
  for (std::size_t i = 0; i < pre_encryption.size(); ++i) {
    EXPECT_NEAR(std::norm(pre_encryption[i] - post_decryption[i]),
                0,
                post_encryption_epsilon);
  }
}

TEST_P(TestPtxtCKKS, canEncryptAndDecryptRealPtxtsWithEa)
{
  helib::buildModChain(context, 100, 2);
  const helib::EncryptedArrayCx& ea = context.ea->getCx();
  helib::SecKey secret_key(context);
  secret_key.GenSecKey();
  const helib::PubKey& public_key(secret_key);

  std::vector<double> data(context.ea->size());
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = (i - 3) / 10.0;
  }
  helib::Ptxt<helib::CKKS> pre_encryption(context, data);
  helib::Ctxt ctxt(public_key);

  public_key.Encrypt(ctxt, pre_encryption);

  helib::Ptxt<helib::CKKS> post_decryption(context);
  ea.decrypt(ctxt, secret_key, post_decryption);
  EXPECT_EQ(pre_encryption.size(), post_decryption.size());
  for (std::size_t i = 0; i < pre_encryption.size(); ++i) {
    EXPECT_NEAR(pre_encryption[i].real(),
                post_decryption[i].real(),
                post_encryption_epsilon);
    EXPECT_NEAR(pre_encryption[i].imag(),
                post_decryption[i].imag(),
                post_encryption_epsilon);
  }
}

TEST_P(TestPtxtCKKS, plusEqualsWithCiphertextWorks)
{
  helib::buildModChain(context, 150, 2);
  helib::SecKey secret_key(context);
  secret_key.GenSecKey();
  const helib::PubKey& public_key(secret_key);

  // Encrypt the augend, addend is plaintext
  std::vector<std::complex<double>> augend_data(context.ea->size());
  std::vector<std::complex<double>> addend_data(context.ea->size());
  for (long i = 0; i < helib::lsize(augend_data); ++i) {
    augend_data[i] = {i / 10.0, -i * i / 63.0};
    addend_data[i] = {-i / 20.0, i * i * i * 2.6};
  }
  helib::Ptxt<helib::CKKS> augend_ptxt(context, augend_data);
  helib::Ptxt<helib::CKKS> addend(context, addend_data);
  helib::Ctxt augend(public_key);
  public_key.Encrypt(augend, augend_ptxt);

  augend += addend;
  augend_ptxt += addend;

  // augend_ptxt and augend should now match
  helib::Ptxt<helib::CKKS> result(context);
  secret_key.Decrypt(result, augend);
  EXPECT_EQ(result.size(), augend_ptxt.size());
  for (std::size_t i = 0; i < result.size(); ++i) {
    EXPECT_NEAR(std::norm(result[i] - augend_ptxt[i]),
                0,
                post_encryption_epsilon);
  }
}

TEST_P(TestPtxtCKKS, addConstantCKKSWithCiphertextWorks)
{
  helib::buildModChain(context, 150, 2);
  helib::SecKey secret_key(context);
  secret_key.GenSecKey();
  const helib::PubKey& public_key(secret_key);

  // Encrypt the augend, addend is plaintext
  std::vector<std::complex<double>> augend_data(context.ea->size());
  std::vector<std::complex<double>> addend_data(context.ea->size());
  for (long i = 0; i < helib::lsize(augend_data); ++i) {
    augend_data[i] = {i / 70.0, -i * 10.5};
    addend_data[i] = {-i / 10.0, i * 0.8};
  }
  helib::Ptxt<helib::CKKS> augend_ptxt(context, augend_data);
  helib::Ptxt<helib::CKKS> addend(context, addend_data);
  helib::Ctxt augend(public_key);
  public_key.Encrypt(augend, augend_ptxt);

  augend.addConstantCKKS(addend);
  augend_ptxt += addend;

  // augend_ptxt and augend should now match
  helib::Ptxt<helib::CKKS> result(context);
  secret_key.Decrypt(result, augend);
  EXPECT_EQ(result.size(), augend_ptxt.size());
  for (std::size_t i = 0; i < result.size(); ++i) {
    EXPECT_NEAR(std::norm(result[i] - augend_ptxt[i]),
                0,
                post_encryption_epsilon);
  }
}

TEST_P(TestPtxtCKKS, minusEqualsWithCiphertextWorks)
{
  helib::buildModChain(context, 150, 2);
  helib::SecKey secret_key(context);
  secret_key.GenSecKey();
  const helib::PubKey& public_key(secret_key);

  // Encrypt the minuend, subtrahend is plaintext
  std::vector<std::complex<double>> minuend_data(context.ea->size());
  std::vector<std::complex<double>> subtrahend_data(context.ea->size());
  for (long i = 0; i < helib::lsize(minuend_data); ++i) {
    minuend_data[i] = {i * i / 30.0, i * i / 4.5};
    subtrahend_data[i] = {(i + 3) / 4.0, -i * i / 1.3};
  }
  helib::Ptxt<helib::CKKS> minuend_ptxt(context, minuend_data);
  helib::Ptxt<helib::CKKS> subtrahend(context, subtrahend_data);
  helib::Ctxt minuend(public_key);
  public_key.Encrypt(minuend, minuend_ptxt);

  minuend -= subtrahend;
  minuend_ptxt -= subtrahend;

  // minuend_ptxt and minuend should now match
  helib::Ptxt<helib::CKKS> result(context);
  secret_key.Decrypt(result, minuend);
  EXPECT_EQ(result.size(), minuend_ptxt.size());
  for (std::size_t i = 0; i < result.size(); ++i) {
    EXPECT_NEAR(std::norm(result[i] - minuend_ptxt[i]),
                0,
                post_encryption_epsilon);
  }
}

TEST_P(TestPtxtCKKS, multByConstantCKKSFromCiphertextWorks)
{
  helib::buildModChain(context, 150, 2);
  helib::SecKey secret_key(context);
  secret_key.GenSecKey();
  const helib::PubKey& public_key(secret_key);

  // Encrypt the multiplier, multiplicand is plaintext
  std::vector<std::complex<double>> multiplier_data(context.ea->size());
  std::vector<std::complex<double>> multiplicand_data(context.ea->size());
  for (long i = 0; i < helib::lsize(multiplier_data); ++i) {
    multiplier_data[i] = {i * 4.5, -i * i / 12.5};
    multiplicand_data[i] = {(i - 2.5) / 3.5, i * 4.2};
  }
  helib::Ptxt<helib::CKKS> multiplier_ptxt(context, multiplier_data);
  helib::Ptxt<helib::CKKS> multiplicand(context, multiplicand_data);

  helib::Ctxt multiplier(public_key);
  public_key.Encrypt(multiplier, multiplier_ptxt);

  multiplier.multByConstantCKKS(multiplicand);
  multiplier_ptxt *= multiplicand;

  // multiplier_ptxt and multiplier should now match
  helib::Ptxt<helib::CKKS> result(context);
  secret_key.Decrypt(result, multiplier);
  EXPECT_EQ(result.size(), multiplier_ptxt.size());
  for (std::size_t i = 0; i < result.size(); ++i) {
    EXPECT_NEAR(std::norm(result[i] - multiplier_ptxt[i]),
                0,
                post_encryption_epsilon);
  }
}

TEST_P(TestPtxtCKKS, timesEqualsFromCiphertextWorks)
{
  helib::buildModChain(context, 150, 2);
  helib::SecKey secret_key(context);
  secret_key.GenSecKey();
  const helib::PubKey& public_key(secret_key);

  // Encrypt the multiplier, multiplicand is plaintext
  std::vector<std::complex<double>> multiplier_data(context.ea->size());
  std::vector<std::complex<double>> multiplicand_data(context.ea->size());
  for (long i = 0; i < helib::lsize(multiplier_data); ++i) {
    multiplier_data[i] = {i * 4.5, -i * i / 3.3};
    multiplicand_data[i] = {(i - 2.5) / 3.5, i * i / 12.4};
  }
  helib::Ptxt<helib::CKKS> multiplier_ptxt(context, multiplier_data);
  helib::Ptxt<helib::CKKS> multiplicand(context, multiplicand_data);

  helib::Ctxt multiplier(public_key);
  public_key.Encrypt(multiplier, multiplier_ptxt);

  multiplier *= multiplicand;
  multiplier_ptxt *= multiplicand;

  // multiplier_ptxt and multiplier should now match
  helib::Ptxt<helib::CKKS> result(context);
  secret_key.Decrypt(result, multiplier);
  EXPECT_EQ(result.size(), multiplier_ptxt.size());
  for (std::size_t i = 0; i < result.size(); ++i) {
    EXPECT_NEAR(std::norm(result[i] - multiplier_ptxt[i]),
                0,
                post_encryption_epsilon);
  }
}

TEST_P(TestPtxtCKKS, plusOperatorWithOtherPtxtWorks)
{
  std::vector<std::complex<double>> augend_data(context.ea->size());
  std::vector<std::complex<double>> addend_data(context.ea->size());
  std::vector<std::complex<double>> expected_sum_data(context.ea->size());
  for (long i = 0; i < helib::lsize(augend_data); ++i) {
    augend_data[i] = {i / 10.0, -i * i / 3.0};
    addend_data[i] = {-i / 20.0, i * i * i * 42.6};
    expected_sum_data[i] = augend_data[i] + addend_data[i];
  }
  helib::Ptxt<helib::CKKS> augend(context, augend_data);
  helib::Ptxt<helib::CKKS> addend(context, addend_data);
  helib::Ptxt<helib::CKKS> sum;

  sum = augend + addend;

  EXPECT_EQ(expected_sum_data.size(), sum.size());
  for (std::size_t i = 0; i < sum.size(); ++i) {
    EXPECT_EQ(expected_sum_data[i], sum[i]);
  }
}

TEST_P(TestPtxtCKKS, minusOperatorWithOtherPtxtWorks)
{
  std::vector<std::complex<double>> minuend_data(context.ea->size());
  std::vector<std::complex<double>> subtrahend_data(context.ea->size());
  std::vector<std::complex<double>> expected_diff_data(context.ea->size());
  for (long i = 0; i < helib::lsize(minuend_data); ++i) {
    minuend_data[i] = {i / 10.0, -i * i / 3.0};
    subtrahend_data[i] = {-i / 20.0, i * i * i * 42.6};
    expected_diff_data[i] = minuend_data[i] - subtrahend_data[i];
  }
  helib::Ptxt<helib::CKKS> minuend(context, minuend_data);
  helib::Ptxt<helib::CKKS> subtrahend(context, subtrahend_data);
  helib::Ptxt<helib::CKKS> diff;

  diff = minuend - subtrahend;

  EXPECT_EQ(expected_diff_data.size(), diff.size());
  for (std::size_t i = 0; i < diff.size(); ++i) {
    EXPECT_EQ(expected_diff_data[i], diff[i]);
  }
}

TEST_P(TestPtxtCKKS, timesOperatorWithOtherPtxtWorks)
{
  std::vector<std::complex<double>> multiplier_data(context.ea->size());
  std::vector<std::complex<double>> multiplicand_data(context.ea->size());
  std::vector<std::complex<double>> expected_product_data(context.ea->size());
  for (long i = 0; i < helib::lsize(multiplier_data); ++i) {
    multiplier_data[i] = {i / 10.0, -i * i / 3.0};
    multiplicand_data[i] = {-i / 20.0, i * i * i * 42.6};
    expected_product_data[i] = multiplier_data[i] * multiplicand_data[i];
  }
  helib::Ptxt<helib::CKKS> multiplier(context, multiplier_data);
  helib::Ptxt<helib::CKKS> multiplicand(context, multiplicand_data);
  helib::Ptxt<helib::CKKS> product;

  product = multiplier * multiplicand;

  EXPECT_EQ(expected_product_data.size(), product.size());
  for (std::size_t i = 0; i < product.size(); ++i) {
    EXPECT_EQ(expected_product_data[i], product[i]);
  }
}

class TestPtxtBGV : public ::testing::TestWithParam<BGVParameters>
{
protected:
  TestPtxtBGV() :
      m(GetParam().m),
      p(GetParam().p),
      r(GetParam().r),
      ppowr(power(p, r)),
      context(m, p, r)
  {}

  static long power(long base, unsigned long exponent)
  {
    long result = base;
    while (--exponent)
      result *= base;
    return result;
  }

  const unsigned long m;
  const unsigned long p;
  const unsigned long r;
  const unsigned long ppowr;

  helib::Context context;
};

TEST_P(TestPtxtBGV, canBeConstructedWithBGVContext)
{
  helib::Ptxt<helib::BGV> ptxt(context);
}

TEST_P(TestPtxtBGV, canBeDefaultConstructed) { helib::Ptxt<helib::BGV> ptxt; }

TEST_P(TestPtxtBGV, canBeCopyConstructed)
{
  helib::Ptxt<helib::BGV> ptxt(context);
  helib::Ptxt<helib::BGV> ptxt2(ptxt);
}

TEST_P(TestPtxtBGV, canBeAssignedFromOtherPtxt)
{
  helib::Ptxt<helib::BGV> ptxt(context);
  helib::Ptxt<helib::BGV> ptxt2 = ptxt;
}

TEST_P(TestPtxtBGV, reportsWhetherItIsValid)
{
  helib::Ptxt<helib::BGV> invalid_ptxt;
  helib::Ptxt<helib::BGV> valid_ptxt(context);
  EXPECT_FALSE(invalid_ptxt.isValid());
  EXPECT_TRUE(valid_ptxt.isValid());
}

TEST_P(TestPtxtBGV, preservesLongDataPassedIntoConstructor)
{
  std::vector<long> data(context.ea->size());
  std::iota(data.begin(), data.end(), 0);
  helib::Ptxt<helib::BGV> ptxt(context, data);
  EXPECT_EQ(ptxt.size(), data.size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    EXPECT_EQ(ptxt[i], data[i]);
  }
}

TEST_P(TestPtxtBGV, preservesCoefficientVectorDataPassedIntoConstructor)
{
  std::vector<std::vector<long>> data(context.ea->size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {1};
  }
  helib::Ptxt<helib::BGV> ptxt(context, data);
  EXPECT_EQ(ptxt.size(), data.size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    EXPECT_EQ(ptxt[i], data[i]);
  }
}

TEST_P(TestPtxtBGV, preservesZzxDataPassedIntoConstructor)
{
  std::vector<NTL::ZZX> data(context.ea->size());
  std::iota(data.begin(), data.end(), 0);
  helib::Ptxt<helib::BGV> ptxt(context, data);
  EXPECT_EQ(ptxt.size(), data.size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    EXPECT_EQ(ptxt[i], data[i]);
  }
}

TEST_P(TestPtxtBGV, writesDataCorrectlyToOstream)
{
  const long p2r = context.slotRing->p2r;
  const long d = context.zMStar.getOrdP();
  helib::PolyMod poly(context.slotRing);
  std::vector<helib::PolyMod> data(context.ea->size(), poly);
  std::stringstream ss;
  ss << "[";
  for (long i = 0; i < helib::lsize(data); ++i) {
    NTL::ZZX input;
    NTL::SetCoeff(input, 0, i % p2r);
    if (d != 1) {
      NTL::SetCoeff(input, 1, (i + 2) % p2r);
    }
    data[i] = input;
    // Serialisation of data[i] (i.e. PolyMod) is tested in `TestPolyMod.cpp`
    ss << data[i] << (i != helib::lsize(data) - 1 ? ", " : "");
  }
  ss << "]";
  helib::Ptxt<helib::BGV> ptxt(context, data);
  std::string expected = ss.str();
  std::ostringstream os;
  os << ptxt;

  EXPECT_EQ(os.str(), expected);
}

TEST_P(TestPtxtBGV, readsDataCorrectlyFromIstream)
{
  helib::PolyMod poly(context.slotRing);
  std::vector<helib::PolyMod> data(context.ea->size(), poly);
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = {i, i + 2};
  }
  helib::Ptxt<helib::BGV> ptxt(context);
  std::stringstream ss;
  ss << "[";
  for (auto it = data.begin(); it != data.end(); it++) {
    ss << *it;
    if (it != data.end() - 1) {
      ss << ", ";
    }
  }
  ss << "]";

  ss >> ptxt;

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], data[i]);
  }
}

TEST_P(TestPtxtBGV, deserializeIsInverseOfSerialize)
{
  helib::PolyMod poly(context.slotRing);
  std::vector<helib::PolyMod> data(context.ea->size(), poly);
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = {i, i + 2};
  }
  helib::Ptxt<helib::BGV> ptxt(context);

  std::stringstream str;
  str << ptxt;

  helib::Ptxt<helib::BGV> deserialized(context);
  str >> deserialized;

  EXPECT_EQ(ptxt, deserialized);
}

TEST_P(TestPtxtBGV, serializeFunctionSerializesCorrectly)
{
  helib::PolyMod poly(context.slotRing);
  std::vector<helib::PolyMod> data(context.ea->size(), poly);
  std::stringstream ptxt_string_stream;
  ptxt_string_stream << "[";
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = 2 * i;
    ptxt_string_stream << "[" << helib::mcMod(2 * i, ppowr) << "]";
    if (i < helib::lsize(data) - 1)
      ptxt_string_stream << ", ";
  }
  ptxt_string_stream << "]";
  helib::Ptxt<helib::BGV> ptxt(context, data);

  std::stringstream ss;
  helib::serialize(ss, ptxt);

  EXPECT_EQ(ss.str(), ptxt_string_stream.str());
}

TEST_P(TestPtxtBGV, deserializeFunctionDeserializesCorrectly)
{
  helib::PolyMod poly(context.slotRing);
  std::vector<helib::PolyMod> data(context.ea->size(), poly);
  std::stringstream ptxt_string_stream;
  ptxt_string_stream << "[";
  for (long i = 0; i < helib::lsize(data); ++i) {
    NTL::ZZX tmp;
    ptxt_string_stream << "[";
    for (long j = 0; j < context.zMStar.getOrdP(); ++j) {
      NTL::SetCoeff(tmp, j, j * j);
      ptxt_string_stream << j * j;
      if (j < context.zMStar.getOrdP() - 1)
        ptxt_string_stream << ",";
    }
    data[i] = tmp;
    ptxt_string_stream << "]";
    if (i < helib::lsize(data) - 1)
      ptxt_string_stream << ", ";
  }
  ptxt_string_stream << "]";
  helib::Ptxt<helib::BGV> ptxt(context, data);

  helib::Ptxt<helib::BGV> deserialized_ptxt(context);
  helib::deserialize(ptxt_string_stream, deserialized_ptxt);

  EXPECT_EQ(ptxt, deserialized_ptxt);
}

TEST_P(TestPtxtBGV, deserializeFunctionThrowsIfMoreElementsThanSlots)
{
  helib::PolyMod poly(context.slotRing);
  std::vector<helib::PolyMod> data(context.ea->size() + 1, poly);
  std::stringstream ptxt_string_stream;
  ptxt_string_stream << "[";
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = {i, i * i};
    ptxt_string_stream << "[" << i << ", " << i * i << "]";
    if (i < helib::lsize(data) - 1)
      ptxt_string_stream << ", ";
  }
  ptxt_string_stream << "]";

  helib::Ptxt<helib::BGV> deserialized_ptxt(context);

  EXPECT_THROW(helib::deserialize(ptxt_string_stream, deserialized_ptxt),
               helib::IOError);
}

TEST_P(TestPtxtBGV, rightShiftOperatorThrowsIfMoreElementsThanSlots)
{
  helib::PolyMod poly(context.slotRing);
  std::vector<helib::PolyMod> data(context.ea->size() + 1, poly);
  std::stringstream ptxt_string_stream;
  ptxt_string_stream << "[";
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = {i, i * i};
    ptxt_string_stream << "[" << i << ", " << i * i << "]";
    if (i < helib::lsize(data) - 1)
      ptxt_string_stream << ", ";
  }
  ptxt_string_stream << "]";

  helib::Ptxt<helib::BGV> deserialized_ptxt(context);
  EXPECT_THROW(ptxt_string_stream >> deserialized_ptxt, helib::IOError);
}

TEST_P(TestPtxtBGV, readsManyPtxtsFromStream)
{
  helib::PolyMod poly(context.slotRing);
  std::vector<helib::PolyMod> data1(context.ea->size(), poly);
  std::vector<helib::PolyMod> data2(context.ea->size(), poly);
  std::vector<helib::PolyMod> data3(context.ea->size(), poly);
  for (long i = 0; i < helib::lsize(data1); ++i) {
    data1[i] = {i, i + 2};
    data2[i] = {2 * i, 2 * (i + 2)};
    data3[i] = {3 * i, 3 * (i + 2)};
  }
  helib::Ptxt<helib::BGV> ptxt1(context, data1);
  helib::Ptxt<helib::BGV> ptxt2(context, data2);
  helib::Ptxt<helib::BGV> ptxt3(context, data3);

  std::stringstream ss;
  ss << ptxt1 << std::endl;
  ss << ptxt2 << std::endl;
  ss << ptxt3 << std::endl;

  helib::Ptxt<helib::BGV> deserialized1(context);
  helib::Ptxt<helib::BGV> deserialized2(context);
  helib::Ptxt<helib::BGV> deserialized3(context);
  ss >> deserialized1;
  ss >> deserialized2;
  ss >> deserialized3;

  EXPECT_EQ(ptxt1, deserialized1);
  EXPECT_EQ(ptxt2, deserialized2);
  EXPECT_EQ(ptxt3, deserialized3);
}

TEST_P(TestPtxtBGV, preservesPolyModDataPassedIntoConstructor)
{
  helib::PolyMod poly(context.slotRing);
  std::vector<helib::PolyMod> data(context.ea->size(), poly);
  helib::Ptxt<helib::BGV> ptxt(context, data);
  EXPECT_EQ(ptxt.size(), data.size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    EXPECT_EQ(ptxt[i], data[i]);
  }
}

TEST_P(TestPtxtBGV, throwsIfp2rAndGDoNotMatchThoseFromContext)
{
  NTL::ZZX G = context.slotRing->G;
  long p = context.slotRing->p;
  long r = context.slotRing->r;
  // Non-matching p^r
  std::shared_ptr<helib::PolyModRing> badPolyModRing1(
      new helib::PolyModRing(p + 1, r, G));
  helib::PolyMod badPolyMod1(badPolyModRing1);
  // Non-matching G
  std::shared_ptr<helib::PolyModRing> badPolyModRing2(
      new helib::PolyModRing(p, r, G + 1));
  helib::PolyMod badPolyMod2(badPolyModRing2);
  // All good
  std::shared_ptr<helib::PolyModRing> goodPolyModRing(
      new helib::PolyModRing(p, r, G));
  helib::PolyMod goodPolyMod(goodPolyModRing);

  std::vector<helib::PolyMod> data(context.ea->size(), goodPolyMod);

  // Make all of them good except 1, make sure it still notices
  data.back() = badPolyMod1;
  EXPECT_THROW(helib::Ptxt<helib::BGV> ptxt(context, data),
               helib::RuntimeError);
  data.back() = badPolyMod2;
  EXPECT_THROW(helib::Ptxt<helib::BGV> ptxt(context, data),
               helib::RuntimeError);

  // Make sure it complains if it's just given 1 bad PolyMod too
  EXPECT_THROW(helib::Ptxt<helib::BGV> ptxt(context, badPolyMod1),
               helib::RuntimeError);
  EXPECT_THROW(helib::Ptxt<helib::BGV> ptxt(context, badPolyMod2),
               helib::RuntimeError);
}

TEST_P(TestPtxtBGV, lsizeReportsCorrectSize)
{
  std::vector<long> data(context.ea->size());
  helib::Ptxt<helib::BGV> ptxt(context, data);
  EXPECT_EQ(ptxt.lsize(), data.size());
}

TEST_P(TestPtxtBGV, sizeReportsCorrectSize)
{
  std::vector<long> data(context.ea->size());
  helib::Ptxt<helib::BGV> ptxt(context, data);
  EXPECT_EQ(ptxt.size(), data.size());
}

TEST_P(TestPtxtBGV, padsWithZerosWhenPassingInSmallDataVector)
{
  std::vector<long> data(context.ea->size() - 1);
  std::iota(data.begin(), data.end(), 0);
  helib::Ptxt<helib::BGV> ptxt(context, data);
  for (std::size_t i = 0; i < data.size(); ++i) {
    EXPECT_EQ(ptxt[i], helib::mcMod(data[i], ppowr));
  }
  for (std::size_t i = data.size(); i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], 0l);
  }
}

TEST_P(TestPtxtBGV, hasSameNumberOfSlotsAsContext)
{
  helib::Ptxt<helib::BGV> ptxt(context);
  EXPECT_EQ(context.ea->size(), ptxt.size());
}

TEST_P(TestPtxtBGV, randomSetsDataRandomly)
{
  helib::Ptxt<helib::BGV> ptxt(context);
  ptxt.random();
  std::vector<helib::Ptxt<helib::BGV>> ptxts(5, ptxt);
  for (auto& p : ptxts)
    p.random();

  bool all_equal = true;
  for (std::size_t i = 0; i < ptxts.size() - 1; ++i)
    if (ptxts[i] != ptxts[i + 1]) {
      all_equal = false;
      break;
    }
  EXPECT_FALSE(all_equal) << "5 random ptxts are all equal - likely that"
                             " random() is not actually randomising!";
}

TEST_P(TestPtxtBGV, runningSumsWorksCorrectly)
{
  std::vector<long> data(context.ea->size());
  std::iota(data.begin(), data.end(), 1);
  std::vector<long> expected_result(data.size());
  for (std::size_t i = 0; i < data.size(); ++i)
    expected_result[i] = ((i + 1) * (i + 2)) / 2;

  helib::Ptxt<helib::BGV> ptxt(context, data);
  ptxt.runningSums();

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], expected_result[i]);
  }
}

TEST_P(TestPtxtBGV, totalSumsWorksCorrectly)
{
  std::vector<long> data(context.ea->size());
  std::iota(data.begin(), data.end(), 1);
  std::vector<long> expected_result(data.size());
  for (std::size_t i = 0; i < data.size(); ++i)
    expected_result[i] = (data.size() * (data.size() + 1)) / 2;

  helib::Ptxt<helib::BGV> ptxt(context, data);
  ptxt.totalSums();

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], expected_result[i]);
  }
}

TEST_P(TestPtxtBGV, incrementalProductWorksCorrectly)
{
  std::vector<long> data(context.ea->size());
  std::iota(data.begin(), data.end(), 1);
  std::vector<long> expected_result(data);
  for (std::size_t i = 1; i < data.size(); ++i)
    expected_result[i] =
        (expected_result[i] * expected_result[i - 1]) % context.slotRing->p2r;

  helib::Ptxt<helib::BGV> ptxt(context, data);
  ptxt.incrementalProduct();

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], expected_result[i]);
  }
}

TEST_P(TestPtxtBGV, totalProductWorksCorrectly)
{
  std::vector<long> data(context.ea->size());
  std::iota(data.begin(), data.end(), 1);
  long product = 1;
  for (std::size_t i = 0; i < data.size(); ++i) {
    product *= data[i];
    product %= context.slotRing->p2r;
  }
  std::vector<long> expected_result(data.size(), product);

  helib::Ptxt<helib::BGV> ptxt(context, data);
  ptxt.totalProduct();

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], expected_result[i]);
  }
}

TEST_P(TestPtxtBGV, innerProductWorksCorrectly)
{
  std::vector<long> data(context.ea->size());
  std::iota(data.begin(), data.end(), 0);
  helib::Ptxt<helib::BGV> ptxt(context, data);
  std::vector<helib::Ptxt<helib::BGV>> first_ptxt_vector(4, ptxt);
  ptxt += ptxt;
  std::vector<helib::Ptxt<helib::BGV>> second_ptxt_vector(4, ptxt);

  helib::Ptxt<helib::BGV> result(context);
  innerProduct(result, first_ptxt_vector, second_ptxt_vector);

  std::vector<long> expected_result(data.size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    expected_result[i] = 4 * (data[i] * (2 * data[i]));
  }

  for (std::size_t i = 0; i < result.size(); ++i) {
    EXPECT_EQ(result[i], expected_result[i]);
  }
}

TEST_P(TestPtxtBGV, mapTo01MapsSlotsCorrectly)
{
  std::vector<long> data(context.ea->size());
  std::iota(data.begin(), data.end(), 0);
  std::vector<long> expected_result(data.size(), 1);
  for (std::size_t i = 0; i < data.size(); ++i)
    if (i % p == 0)
      expected_result[i] = 0;

  // Should exist as a free function and a member function
  helib::Ptxt<helib::BGV> ptxt(context, data);
  helib::Ptxt<helib::BGV> ptxt2(context, data);
  ptxt.mapTo01();
  mapTo01(*(context.ea), ptxt2);

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], expected_result[i]);
    EXPECT_EQ(ptxt2[i], expected_result[i]);
  }
}

TEST(TestPtxtBGV, automorphWorksCorrectly)
{
  std::vector<long> gens = {11, 2};
  std::vector<long> ords = {6, 2};
  const helib::Context context(45, 19, 1, gens, ords);
  std::vector<NTL::ZZX> data(context.ea->size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    NTL::SetX(data[i]);
    (data[i] += 1) *= i;
  }

  helib::Ptxt<helib::BGV> ptxt(context, data);
  helib::Ptxt<helib::BGV> expected_result(ptxt);
  expected_result[0] = {13, 10};
  expected_result[1] = {0, 0};
  expected_result[2] = {18, 8};
  expected_result[3] = {2, 2};
  expected_result[4] = {12, 18};
  expected_result[5] = {4, 4};
  expected_result[6] = {17, 16};
  expected_result[7] = {6, 6};
  expected_result[8] = {3, 14};
  expected_result[9] = {8, 8};
  expected_result[10] = {8, 12};
  expected_result[11] = {10, 10};

  ptxt.automorph(2);

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], expected_result[i]);
  }
}

TEST_P(TestPtxtBGV, frobeniusAutomorphWithConstantsWorksCorrectly)
{
  std::vector<long> data(context.ea->size());
  std::iota(data.begin(), data.end(), 0);
  std::vector<long> expected_result(data);
  helib::Ptxt<helib::BGV> ptxt(context, data);
  for (long i = 0; i <= context.zMStar.getOrdP(); ++i) {
    auto ptxtUnderTest = ptxt;
    ptxtUnderTest.frobeniusAutomorph(i);
    for (std::size_t j = 0; j < ptxtUnderTest.size(); ++j) {
      ASSERT_EQ(ptxtUnderTest[j], expected_result[j])
          << "i = " << i << " j = " << j << std::endl;
    }
  }
}

TEST(TestPtxtBGV, frobeniusAutomorphWithPolynomialsWorksCorrectly)
{
  const helib::Context context(45, 19, 1);
  std::vector<NTL::ZZX> data(context.ea->size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    NTL::SetX(data[i]);
    (data[i] += 1) *= i;
  }

  helib::Ptxt<helib::BGV> ptxt(context, data);
  helib::Ptxt<helib::BGV> expected_result(ptxt);
  expected_result[0] = {0, 0};
  expected_result[1] = {12, 18};
  expected_result[2] = {5, 17};
  expected_result[3] = {17, 16};
  expected_result[4] = {10, 15};
  expected_result[5] = {3, 14};
  expected_result[6] = {15, 13};
  expected_result[7] = {8, 12};
  expected_result[8] = {1, 11};
  expected_result[9] = {13, 10};
  expected_result[10] = {6, 9};
  expected_result[11] = {18, 8};

  ptxt.frobeniusAutomorph(1);

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], expected_result[i]);
  }
}

TEST_P(TestPtxtBGV, timesEqualsOtherPlaintextWorks)
{
  std::vector<long> product_data(context.ea->size(), 3);
  std::vector<long> multiplier_data(context.ea->size());
  std::iota(multiplier_data.begin(), multiplier_data.end(), 0);

  std::vector<long> expected_result(product_data);
  for (std::size_t i = 0; i < product_data.size(); ++i) {
    expected_result[i] = expected_result[i] * multiplier_data[i];
  }

  helib::Ptxt<helib::BGV> product(context, product_data);
  helib::Ptxt<helib::BGV> multiplier(context, multiplier_data);

  product *= multiplier;

  for (std::size_t i = 0; i < product.size(); ++i) {
    EXPECT_EQ(product[i], expected_result[i]);
  }
}

TEST_P(TestPtxtBGV, minusEqualsOtherPlaintextWorks)
{
  std::vector<long> difference_data(context.ea->size(), 1);
  std::vector<long> subtrahend_data(context.ea->size());
  std::iota(subtrahend_data.begin(), subtrahend_data.end(), 0);

  std::vector<long> expected_result(difference_data);
  for (std::size_t i = 0; i < subtrahend_data.size(); ++i) {
    expected_result[i] = expected_result[i] - subtrahend_data[i];
  }

  helib::Ptxt<helib::BGV> difference(context, difference_data);
  helib::Ptxt<helib::BGV> subtrahend(context, subtrahend_data);

  difference -= subtrahend;

  for (std::size_t i = 0; i < difference.size(); ++i) {
    EXPECT_EQ(difference[i], expected_result[i]);
  }
}

TEST_P(TestPtxtBGV, minusEqualsScalarWorks)
{
  std::vector<long> data(context.ea->size());
  std::iota(data.begin(), data.end(), 0);

  const long scalar = 3;

  std::vector<long> expected_result(data);
  for (auto& num : expected_result)
    num = num - scalar;

  helib::Ptxt<helib::BGV> ptxt(context, data);
  ptxt -= scalar;

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], expected_result[i]);
  }
}

TEST_P(TestPtxtBGV, plusEqualsOtherPlaintextWorks)
{
  std::vector<long> augend_data(context.ea->size());
  std::iota(augend_data.begin(), augend_data.end(), 0);
  std::vector<long> addend_data(context.ea->size());
  for (long i = 0; i < helib::lsize(addend_data); ++i)
    addend_data[i] = helib::mcMod(2 * i + 1, p);
  std::vector<long> expected_result(context.ea->size());
  for (long i = 0; i < helib::lsize(expected_result); ++i)
    expected_result[i] = augend_data[i] + addend_data[i];

  helib::Ptxt<helib::BGV> sum(context, augend_data);
  helib::Ptxt<helib::BGV> addend(context, addend_data);
  sum += addend;

  for (std::size_t i = 0; i < sum.size(); ++i) {
    EXPECT_EQ(sum[i], expected_result[i]);
  }
}

TEST_P(TestPtxtBGV, plusEqualsScalarWorks)
{
  std::vector<long> data(context.ea->size());
  std::iota(data.begin(), data.end(), 0);

  const long scalar = 3;

  std::vector<long> expected_result(data);
  for (auto& num : expected_result)
    num = num + scalar;

  helib::Ptxt<helib::BGV> ptxt(context, data);
  ptxt += scalar;

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], expected_result[i]);
  }
}

TEST_P(TestPtxtBGV, timesEqualsScalarWorks)
{
  std::vector<long> data(context.ea->size());
  std::iota(data.begin(), data.end(), 0);

  const long scalar = 3;

  std::vector<long> expected_result(data);
  for (auto& num : expected_result)
    num = num * scalar;

  helib::Ptxt<helib::BGV> ptxt(context, data);
  ptxt *= scalar;

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], expected_result[i]);
  }
}

TEST_P(TestPtxtBGV, equalityWithOtherPlaintextWorks)
{
  std::vector<long> data(context.ea->size());
  std::iota(data.begin(), data.end(), 0);

  helib::Ptxt<helib::BGV> ptxt1(context, data);
  helib::Ptxt<helib::BGV> ptxt2(context, data);
  EXPECT_TRUE(ptxt1 == ptxt2);
}

TEST_P(TestPtxtBGV, notEqualsOperatorWithOtherPlaintextWorks)
{
  std::vector<long> data1(context.ea->size());
  std::iota(data1.begin(), data1.end(), 0);
  std::vector<long> data2(context.ea->size());
  std::iota(data2.begin(), data2.end(), 1);

  helib::Ptxt<helib::BGV> ptxt1(context, data1);
  helib::Ptxt<helib::BGV> ptxt2(context, data2);
  EXPECT_TRUE(ptxt1 != ptxt2);
}

TEST_P(TestPtxtBGV, negateNegatesCorrectly)
{
  std::vector<long> data(context.ea->size());
  std::iota(data.begin(), data.end(), 0);

  std::vector<long> expected_result(data);
  for (auto& num : expected_result)
    num = -num;

  helib::Ptxt<helib::BGV> ptxt(context, data);
  ptxt.negate();

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], expected_result[i]);
  }
}

TEST_P(TestPtxtBGV, addConstantWorksCorrectly)
{
  NTL::ZZX input;
  NTL::SetCoeff(input, 0, 2);
  NTL::SetCoeff(input, 1, 1);

  helib::PolyMod poly(input, context.slotRing);
  std::vector<helib::PolyMod> data(context.ea->size());
  for (std::size_t i = 0; i < data.size(); ++i)
    data[i] = poly + i;

  std::vector<helib::PolyMod> expected_result(data);
  for (auto& num : expected_result)
    (num += input) += 3L;

  helib::Ptxt<helib::BGV> ptxt(context, data);
  ptxt.addConstant(input).addConstant(3l);

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], expected_result[i]);
  }
}

TEST_P(TestPtxtBGV, multiplyByMultipliesCorrectly)
{
  std::vector<long> product_data(context.ea->size(), 3);
  std::vector<long> multiplier_data(context.ea->size());
  std::iota(multiplier_data.begin(), multiplier_data.end(), 0);

  std::vector<long> expected_result(product_data);
  for (std::size_t i = 0; i < product_data.size(); ++i) {
    expected_result[i] = expected_result[i] * multiplier_data[i];
  }

  helib::Ptxt<helib::BGV> product(context, product_data);
  helib::Ptxt<helib::BGV> multiplier(context, multiplier_data);

  product.multiplyBy(multiplier);

  for (std::size_t i = 0; i < product.size(); ++i) {
    EXPECT_EQ(product[i], expected_result[i]);
  }
}

TEST_P(TestPtxtBGV, multiplyBy2MultipliesCorrectly)
{
  std::vector<long> product_data(context.ea->size(), 3);
  std::vector<long> multiplier_data1(context.ea->size());
  std::vector<long> multiplier_data2(context.ea->size());
  std::iota(multiplier_data1.begin(), multiplier_data1.end(), 0);
  std::iota(multiplier_data2.begin(), multiplier_data2.end(), 0);

  std::vector<long> expected_result(product_data);
  for (std::size_t i = 0; i < product_data.size(); ++i) {
    expected_result[i] =
        expected_result[i] * multiplier_data1[i] * multiplier_data2[i];
  }

  helib::Ptxt<helib::BGV> product(context, product_data);
  helib::Ptxt<helib::BGV> multiplier1(context, multiplier_data1);
  helib::Ptxt<helib::BGV> multiplier2(context, multiplier_data2);

  product.multiplyBy2(multiplier1, multiplier2);

  for (std::size_t i = 0; i < product.size(); ++i) {
    EXPECT_EQ(product[i], expected_result[i]);
  }
}

TEST_P(TestPtxtBGV, squareSquaresCorrectly)
{
  std::vector<long> data(context.ea->size());
  std::iota(data.begin(), data.end(), 0);
  std::vector<long> expected_result(data);
  for (auto& num : expected_result)
    num = num * num;
  helib::Ptxt<helib::BGV> ptxt(context, data);
  ptxt.square();
  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], expected_result[i]);
  }
}

TEST_P(TestPtxtBGV, cubeCubesCorrectly)
{
  std::vector<long> data(context.ea->size());
  std::iota(data.begin(), data.end(), 0);
  std::vector<long> expected_result(data);
  for (auto& num : expected_result)
    num = num * num * num;
  helib::Ptxt<helib::BGV> ptxt(context, data);
  ptxt.cube();
  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], expected_result[i]);
  }
}

TEST_P(TestPtxtBGV, powerCorrectlyRaisesToPowers)
{
  std::vector<long> data(context.ea->size());
  std::iota(data.begin(),
            data.end(),
            -(static_cast<long>(context.ea->size()) / 2));
  std::vector<long> exponents{1, 3, 4, 5, 300};

  const auto naive_powermod =
      [](long base, unsigned long exponent, unsigned long mod) {
        if (exponent == 0)
          return 1l;

        long result = base;
        while (--exponent)
          result = helib::mcMod(result * base, mod);
        return result;
      };

  for (const auto& exponent : exponents) {
    std::vector<long> expected_result(data);
    for (auto& num : expected_result)
      num = naive_powermod(num, exponent, ppowr);
    helib::Ptxt<helib::BGV> ptxt(context, data);
    ptxt.power(exponent);
    for (std::size_t i = 0; i < ptxt.size(); ++i) {
      EXPECT_EQ(ptxt[i], expected_result[i]);
    }
  }

  // Make sure raising to 0 throws
  helib::Ptxt<helib::CKKS> ptxt(context, data);
  EXPECT_THROW(ptxt.power(0l), helib::InvalidArgument);
}

TEST_P(TestPtxtBGV, shiftShiftsRightCorrectly)
{
  std::vector<long> data(context.ea->size());
  std::vector<long> right_shifted_data(context.ea->size());
  const auto non_neg_mod = [](int x, int mod) {
    return ((x % mod) + mod) % mod;
  };
  for (std::size_t i = 0; i < data.size(); ++i) {
    if (i > 3) {
      right_shifted_data[i] = non_neg_mod(i - 3, data.size());
    }
    data[i] = i;
  }
  helib::Ptxt<helib::BGV> ptxt(context, data);

  ptxt.shift(3);
  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], right_shifted_data[i]);
  }
}

TEST_P(TestPtxtBGV, shiftShiftsLeftCorrectly)
{
  std::vector<long> data(context.ea->size());
  std::vector<long> left_shifted_data(context.ea->size());
  const auto non_neg_mod = [](int x, int mod) {
    return ((x % mod) + mod) % mod;
  };
  for (std::size_t i = 0; i < data.size(); ++i) {
    if (i < data.size() - 3 && data.size() > 3) {
      left_shifted_data[i] = non_neg_mod(i + 3, data.size());
    }
    data[i] = i;
  }
  helib::Ptxt<helib::BGV> ptxt(context, data);

  ptxt.shift(-3);
  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], left_shifted_data[i]);
  }
}

TEST(TestPtxtBGV, shift1DShiftsRightCorrectly)
{
  long amount = 1;
  const helib::Context context(45, 19, 1);
  std::vector<long> data(context.ea->size());
  std::vector<long> right_shifted_data(context.ea->size());
  const auto shift_first_dim = [](long amount, std::vector<long>& data) {
    std::vector<long> new_data(data.size(), 0l);
    for (long i = 0; i < helib::lsize(data); ++i)
      if (i + 2 * amount < 12 && i + 2 * amount >= 0)
        new_data[i + 2 * amount] = data[i];
    data = std::move(new_data);
  };
  const auto shift_second_dim = [](long amount, std::vector<long>& data) {
    std::vector<long> new_data(data.size(), 0l);
    for (long i = 0; i < helib::lsize(data); ++i)
      switch (amount) {
      case 1l:
        if (i < helib::lsize(data) - 1)
          new_data[i + 1] = i & 1 ? 0 : data[i];
        break;
      case 0l:
        new_data[i] = data[i];
        break;
      case -1l:
        if (i > 0)
          new_data[i - 1] = i & 1 ? data[i] : 0;
        break;
      default:
        new_data[i] = 0;
        break;
      }
    data = std::move(new_data);
  };
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = i;
  }
  {
    helib::Ptxt<helib::BGV> ptxt(context, data);

    right_shifted_data = data;
    shift_first_dim(amount, right_shifted_data);
    ptxt.shift1D(0, amount);

    for (std::size_t i = 0; i < ptxt.size(); ++i) {
      EXPECT_EQ(ptxt[i], right_shifted_data[i]);
    }
  }
  {
    helib::Ptxt<helib::BGV> ptxt(context, data);

    right_shifted_data = data;
    shift_second_dim(amount, right_shifted_data);
    ptxt.shift1D(1, amount);

    for (std::size_t i = 0; i < ptxt.size(); ++i) {
      EXPECT_EQ(ptxt[i], right_shifted_data[i]);
    }
  }
}

TEST(TestPtxtBGV, shift1DShiftsLeftCorrectly)
{
  long amount = -1;
  const helib::Context context(45, 19, 1);
  std::vector<long> data(context.ea->size());
  std::vector<long> right_shifted_data(context.ea->size());
  const auto shift_first_dim = [](long amount, std::vector<long>& data) {
    std::vector<long> new_data(data.size(), 0l);
    for (long i = 0; i < helib::lsize(data); ++i)
      if (i + 2 * amount < 12 && i + 2 * amount >= 0)
        new_data[i + 2 * amount] = data[i];
    data = std::move(new_data);
  };
  const auto shift_second_dim = [](long amount, std::vector<long>& data) {
    std::vector<long> new_data(data.size(), 0l);
    for (long i = 0; i < helib::lsize(data); ++i)
      switch (amount) {
      case 1l:
        if (i < helib::lsize(data) - 1)
          new_data[i + 1] = i & 1 ? 0 : data[i];
        break;
      case 0l:
        new_data[i] = data[i];
        break;
      case -1l:
        if (i > 0)
          new_data[i - 1] = i & 1 ? data[i] : 0;
        break;
      default:
        new_data[i] = 0;
        break;
      }
    data = std::move(new_data);
  };
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = i;
  }
  {
    helib::Ptxt<helib::BGV> ptxt(context, data);

    right_shifted_data = data;
    shift_first_dim(amount, right_shifted_data);
    ptxt.shift1D(0, amount);

    for (std::size_t i = 0; i < ptxt.size(); ++i) {
      EXPECT_EQ(ptxt[i], right_shifted_data[i]);
    }
  }
  {
    helib::Ptxt<helib::BGV> ptxt(context, data);

    right_shifted_data = data;
    shift_second_dim(amount, right_shifted_data);
    ptxt.shift1D(1, amount);

    for (std::size_t i = 0; i < ptxt.size(); ++i) {
      EXPECT_EQ(ptxt[i], right_shifted_data[i]);
    }
  }
}

TEST(TestPtxtBGV, rotate1DRotatesCorrectly)
{
  long amount = 1;
  const helib::Context context(45, 19, 1);
  std::vector<long> data(context.ea->size());
  std::vector<long> left_rotated_data(context.ea->size());
  const auto rotate_first_dim = [](long amount, std::vector<long>& data) {
    amount = helib::mcMod(amount, 12);
    std::vector<long> new_data(data);
    for (long i = 0; i < helib::lsize(data); ++i)
      new_data[(i + 2 * amount) % 12] = data[i];
    data = std::move(new_data);
  };
  const auto rotate_second_dim = [](long amount, std::vector<long>& data) {
    std::vector<long> new_data(data);
    for (long i = 0; i < helib::lsize(data); ++i)
      if (amount % 2)
        new_data[i + (i & 1 ? -1 : 1)] = data[i];
    data = std::move(new_data);
  };
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = i;
  }
  helib::Ptxt<helib::BGV> ptxt(context, data);

  // Rotate in first dimension (Good Dimension)
  rotate_first_dim(-amount, data);
  ptxt.rotate1D(0, -amount);

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], data[i]);
  }

  rotate_first_dim(amount, data);
  ptxt.rotate1D(0, amount);
  // Rotating back and forth gives the original data back
  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], data[i]);
  }

  // Rotate in second dimension (Bad Dimension)
  rotate_second_dim(-amount, data);
  ptxt.rotate1D(1, -amount);

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], data[i]);
  }

  rotate_second_dim(amount, data);
  ptxt.rotate1D(1, amount);
  // Rotating back and forth gives the original data back
  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], data[i]);
  }
}

TEST_P(TestPtxtBGV, rotateRotatesCorrectly)
{
  std::vector<long> data(context.ea->size());
  std::vector<long> left_rotated_data(context.ea->size());
  const auto non_neg_mod = [](int x, int mod) {
    return ((x % mod) + mod) % mod;
  };
  for (int i = 0; i < helib::lsize(data); ++i) {
    data[i] = non_neg_mod(i - 3, data.size());
    left_rotated_data[i] = i;
  }
  helib::Ptxt<helib::BGV> ptxt(context, data);

  ptxt.rotate(-3);
  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], left_rotated_data[i]);
  }
  ptxt.rotate(3);
  // Rotating back and forth gives the original data back
  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], data[i]);
  }
}

TEST_P(TestPtxtBGV, replicateReplicatesCorrectly)
{
  std::vector<long> data(context.ea->size());
  std::iota(data.begin(), data.end(), 0);
  helib::Ptxt<helib::BGV> ptxt(context, data);
  helib::replicate(*context.ea, ptxt, data.size() - 1);
  std::vector<long> replicated_data(context.ea->size(), data[data.size() - 1]);
  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], replicated_data[i]);
  }
}

TEST_P(TestPtxtBGV, replicateAllWorksCorrectly)
{
  std::vector<long> data(context.ea->size());
  std::iota(data.begin(), data.end(), 0);
  helib::Ptxt<helib::BGV> ptxt(context, data);
  std::vector<helib::Ptxt<helib::BGV>> replicated_ptxts = ptxt.replicateAll();
  for (long i = 0; i < helib::lsize(data); ++i) {
    for (const auto& slot : replicated_ptxts[i].getSlotRepr()) {
      EXPECT_EQ(slot, data[i]);
    }
  }
}

TEST_P(TestPtxtBGV, clearZeroesAllSlots)
{
  std::vector<long> data(context.ea->size());
  std::iota(data.begin(), data.end(), 0);
  std::vector<long> expected_result(context.ea->size(), 0);
  helib::Ptxt<helib::BGV> ptxt(context, data);
  ptxt.clear();
  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], expected_result[i]);
  }
}

TEST_P(TestPtxtBGV, defaultConstructedPtxtThrowsWhenOperatedOn)
{
  helib::Ptxt<helib::BGV> ptxt1;
  helib::Ptxt<helib::BGV> ptxt2;
  helib::PolyMod poly;
  EXPECT_THROW(ptxt1.getSlotRepr(), helib::RuntimeError);
  EXPECT_THROW(ptxt1.setData(poly), helib::RuntimeError);
  EXPECT_THROW(ptxt1[0], helib::RuntimeError);
  EXPECT_THROW(ptxt1 *= ptxt2, helib::RuntimeError);
  EXPECT_THROW(ptxt1 += ptxt2, helib::RuntimeError);
  EXPECT_THROW(ptxt1 -= ptxt2, helib::RuntimeError);
  EXPECT_THROW(ptxt1 += 1l, helib::RuntimeError);
  EXPECT_THROW(ptxt1 *= 3l, helib::RuntimeError);
  EXPECT_THROW(ptxt1.negate(), helib::RuntimeError);
  EXPECT_THROW(ptxt1.multiplyBy(ptxt2), helib::RuntimeError);
  EXPECT_THROW(ptxt1.multiplyBy2(ptxt1, ptxt2), helib::RuntimeError);
  EXPECT_THROW(ptxt1.square(), helib::RuntimeError);
  EXPECT_THROW(ptxt1.cube(), helib::RuntimeError);
  EXPECT_THROW(ptxt1.power(4l), helib::RuntimeError);
  EXPECT_THROW(ptxt1.size(), helib::RuntimeError);
  EXPECT_THROW(ptxt1.rotate(1), helib::RuntimeError);
  EXPECT_THROW(ptxt1.rotate1D(0, 1), helib::RuntimeError);
  EXPECT_THROW(ptxt1.shift(1), helib::RuntimeError);
  EXPECT_THROW(ptxt1.lsize(), helib::LogicError);
}

TEST_P(TestPtxtBGV, defaultConstructedContextCannotBeRightOperand)
{
  std::vector<long> data(context.ea->size());
  std::iota(data.begin(), data.end(), 0);
  helib::Ptxt<helib::BGV> valid_ptxt(context, data);
  helib::Ptxt<helib::BGV> invalid_ptxt;

  EXPECT_THROW(valid_ptxt *= invalid_ptxt, helib::RuntimeError);
  EXPECT_THROW(valid_ptxt += invalid_ptxt, helib::RuntimeError);
  EXPECT_THROW(valid_ptxt -= invalid_ptxt, helib::RuntimeError);
  EXPECT_THROW(valid_ptxt.multiplyBy(invalid_ptxt), helib::RuntimeError);
  EXPECT_THROW(valid_ptxt.multiplyBy2(invalid_ptxt, valid_ptxt),
               helib::RuntimeError);
  EXPECT_THROW(valid_ptxt.multiplyBy2(valid_ptxt, invalid_ptxt),
               helib::RuntimeError);
}

TEST_P(TestPtxtBGV, cannotOperateBetweenPtxtsWithDifferentContexts)
{
  helib::Context different_context = helib::Context(m, p, 2 * r);
  std::vector<long> data(context.ea->size(), 1);
  helib::Ptxt<helib::BGV> ptxt1(context, data);
  helib::Ptxt<helib::BGV> ptxt2(different_context, data);
  EXPECT_THROW(ptxt1 *= ptxt2, helib::LogicError);
  EXPECT_THROW(ptxt1 -= ptxt2, helib::LogicError);
  EXPECT_THROW(ptxt1 += ptxt2, helib::LogicError);
  EXPECT_THROW(ptxt1.multiplyBy(ptxt2), helib::LogicError);
  EXPECT_THROW(ptxt1.multiplyBy2(ptxt1, ptxt2), helib::LogicError);
}

TEST_P(TestPtxtBGV, preservesDataPassedAsZZX)
{
  // Put in x + 1 and make sure we get x + 1 out
  NTL::ZZX input_polynomial;
  SetCoeff(input_polynomial, 0, 1);
  SetCoeff(input_polynomial, 1, 1);

  helib::Ptxt<helib::BGV> ptxt(context, input_polynomial);

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], input_polynomial);
  }
}

TEST_P(TestPtxtBGV, setDataWorksWithZZXSameOrderAsPhiMX)
{
  NTL::ZZX phi_mx;
  switch (context.alMod.getTag()) {
  case helib::PA_GF2_tag:
    phi_mx = NTL::conv<NTL::ZZX>(
        context.alMod.getDerived(helib::PA_GF2()).getPhimXMod());
    break;
  case helib::PA_zz_p_tag:
    helib::convert(phi_mx,
                   context.alMod.getDerived(helib::PA_zz_p()).getPhimXMod());
    break;
  case helib::PA_cx_tag:
    // CKKS: do nothing
    break;
  default:
    throw helib::LogicError("No valid tag found in EncryptedArray");
  }
  // Put phi_mx + 1 as data
  NTL::ZZX input_polynomial(phi_mx + 1);

  helib::Ptxt<helib::BGV> ptxt(context, input_polynomial);

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], 1l);
  }
}

TEST_P(TestPtxtBGV, decodeSetDataWorks)
{
  // Put in x + 1 and make sure we get x + 1 out
  NTL::ZZX input_polynomial;
  SetCoeff(input_polynomial, 0, 1);
  SetCoeff(input_polynomial, 1, 1);

  std::vector<NTL::ZZX> test_decoded;
  context.ea->decode(test_decoded, input_polynomial);

  helib::Ptxt<helib::BGV> ptxt(context);

  ptxt.decodeSetData(input_polynomial);

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], test_decoded[i]);
  }
}

TEST_P(TestPtxtBGV, decodeSetDataWorksWithZZXSameOrderAsPhiMX)
{
  NTL::ZZX phi_mx;
  switch (context.alMod.getTag()) {
  case helib::PA_GF2_tag:
    phi_mx = NTL::conv<NTL::ZZX>(
        context.alMod.getDerived(helib::PA_GF2()).getPhimXMod());
    break;
  case helib::PA_zz_p_tag:
    helib::convert(phi_mx,
                   context.alMod.getDerived(helib::PA_zz_p()).getPhimXMod());
    break;
  case helib::PA_cx_tag:
    // CKKS: do nothing
    break;
  default:
    throw helib::LogicError("No valid tag found in EncryptedArray");
  }
  // Put phi_mx + 1 as data
  NTL::ZZX input_polynomial(phi_mx + 1);

  helib::Ptxt<helib::BGV> ptxt(context, input_polynomial);

  std::vector<NTL::ZZX> test_decoded;
  context.ea->decode(test_decoded, input_polynomial);

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], test_decoded[i]);
  }
}

TEST_P(TestPtxtBGV, canEncryptAndDecryptPtxts)
{
  helib::buildModChain(context, 30, 2);
  helib::SecKey secret_key(context);
  secret_key.GenSecKey();
  const helib::PubKey& public_key(secret_key);

  std::vector<long> data(context.ea->size());
  std::iota(data.begin(), data.end(), 0);
  helib::Ptxt<helib::BGV> pre_encryption(context, data);
  helib::Ctxt ctxt(public_key);
  public_key.Encrypt(ctxt, pre_encryption);
  helib::Ptxt<helib::BGV> post_decryption(context);
  secret_key.Decrypt(post_decryption, ctxt);
  for (std::size_t i = 0; i < pre_encryption.size(); ++i) {
    EXPECT_EQ(pre_encryption[i], post_decryption[i]);
  }
}

TEST_P(TestPtxtBGV, plusEqualsWithCiphertextWorks)
{
  helib::buildModChain(context, 30, 2);
  helib::SecKey secret_key(context);
  secret_key.GenSecKey();
  const helib::PubKey& public_key(secret_key);

  // Encrypt the augend, addend is plaintext
  std::vector<long> augend_data(context.ea->size());
  std::vector<long> addend_data(context.ea->size());
  std::iota(augend_data.begin(), augend_data.end(), 0);
  std::iota(addend_data.begin(), addend_data.end(), 7);
  helib::Ptxt<helib::BGV> augend_ptxt(context, augend_data);
  helib::Ptxt<helib::BGV> addend(context, addend_data);
  helib::Ctxt augend(public_key);
  public_key.Encrypt(augend, augend_ptxt);

  augend += addend;
  augend_ptxt += addend;

  // augend_ptxt and augend should now match
  helib::Ptxt<helib::BGV> result(context);
  secret_key.Decrypt(result, augend);
  for (std::size_t i = 0; i < result.size(); ++i) {
    EXPECT_EQ(result[i], augend_ptxt[i]);
  }
}

TEST_P(TestPtxtBGV, addConstantFromCiphertextWorks)
{
  helib::buildModChain(context, 30, 2);
  helib::SecKey secret_key(context);
  secret_key.GenSecKey();
  const helib::PubKey& public_key(secret_key);

  // Encrypt the augend, addend is plaintext
  std::vector<long> augend_data(context.ea->size());
  std::vector<long> addend_data(context.ea->size());
  std::iota(augend_data.begin(), augend_data.end(), 0);
  std::iota(addend_data.begin(), addend_data.end(), 7);
  helib::Ptxt<helib::BGV> augend_ptxt(context, augend_data);
  helib::Ptxt<helib::BGV> addend(context, addend_data);
  helib::Ctxt augend(public_key);
  public_key.Encrypt(augend, augend_ptxt);

  augend.addConstant(addend);
  augend_ptxt += addend;

  // augend_ptxt and augend should now match
  helib::Ptxt<helib::BGV> result(context);
  secret_key.Decrypt(result, augend);
  for (std::size_t i = 0; i < result.size(); ++i) {
    EXPECT_EQ(result[i], augend_ptxt[i]);
  }
}

TEST_P(TestPtxtBGV, minusEqualsWithCiphertextWorks)
{
  helib::buildModChain(context, 30, 2);
  helib::SecKey secret_key(context);
  secret_key.GenSecKey();
  const helib::PubKey& public_key(secret_key);

  // Encrypt the minuend, subtrahend is plaintext
  std::vector<long> minuend_data(context.ea->size());
  std::vector<long> subtrahend_data(context.ea->size());
  std::iota(minuend_data.begin(), minuend_data.end(), 0);
  std::iota(subtrahend_data.begin(), subtrahend_data.end(), 7);
  helib::Ptxt<helib::BGV> minuend_ptxt(context, minuend_data);
  helib::Ptxt<helib::BGV> subtrahend(context, subtrahend_data);
  helib::Ctxt minuend(public_key);
  public_key.Encrypt(minuend, minuend_ptxt);

  minuend -= subtrahend;
  minuend_ptxt -= subtrahend;

  // minuend_ptxt and minuend should now match
  helib::Ptxt<helib::BGV> result(context);
  secret_key.Decrypt(result, minuend);
  for (std::size_t i = 0; i < result.size(); ++i) {
    EXPECT_EQ(result[i], minuend_ptxt[i]);
  }
}

TEST_P(TestPtxtBGV, timesEqualsWithCiphertextWorks)
{
  helib::buildModChain(context, 30, 2);
  helib::SecKey secret_key(context);
  secret_key.GenSecKey();
  const helib::PubKey& public_key(secret_key);

  // Encrypt the multiplier, multiplicand is plaintext
  std::vector<long> multiplier_data(context.ea->size());
  std::vector<long> multiplicand_data(context.ea->size());
  std::iota(multiplier_data.begin(), multiplier_data.end(), 0);
  std::iota(multiplicand_data.begin(), multiplicand_data.end(), 7);
  helib::Ptxt<helib::BGV> multiplier_ptxt(context, multiplier_data);
  helib::Ptxt<helib::BGV> multiplicand(context, multiplicand_data);
  helib::Ctxt multiplier(public_key);
  public_key.Encrypt(multiplier, multiplier_ptxt);

  multiplier *= multiplicand;
  multiplier_ptxt *= multiplicand;

  // multiplier_ptxt and multiplier should now match
  helib::Ptxt<helib::BGV> result(context);
  secret_key.Decrypt(result, multiplier);
  for (std::size_t i = 0; i < result.size(); ++i) {
    EXPECT_EQ(result[i], multiplier_ptxt[i]);
  }
}

TEST_P(TestPtxtBGV, multByConstantFromCiphertextWorks)
{
  helib::buildModChain(context, 30, 2);
  helib::SecKey secret_key(context);
  secret_key.GenSecKey();
  const helib::PubKey& public_key(secret_key);

  // Encrypt the multiplier, multiplicand is plaintext
  std::vector<long> multiplier_data(context.ea->size());
  std::vector<long> multiplicand_data(context.ea->size());
  std::iota(multiplier_data.begin(), multiplier_data.end(), 0);
  std::iota(multiplicand_data.begin(), multiplicand_data.end(), 7);
  helib::Ptxt<helib::BGV> multiplier_ptxt(context, multiplier_data);
  helib::Ptxt<helib::BGV> multiplicand(context, multiplicand_data);
  helib::Ctxt multiplier(public_key);
  public_key.Encrypt(multiplier, multiplier_ptxt);

  multiplier.multByConstant(multiplicand);
  multiplier_ptxt *= multiplicand;

  // multiplier_ptxt and multiplier should now match
  helib::Ptxt<helib::BGV> result(context);
  secret_key.Decrypt(result, multiplier);
  for (std::size_t i = 0; i < result.size(); ++i) {
    EXPECT_EQ(result[i], multiplier_ptxt[i]);
  }
}

// Useful for testing non-power of 2 for CKKS
// INSTANTIATE_TEST_SUITE_P(various_parameters, TestPtxtCKKS,
//        ::testing::Values( 17, 168, 126, 78, 33, 50, 64)
// );

INSTANTIATE_TEST_SUITE_P(
    various_Parameters,
    TestPtxtCKKS,
    ::testing::Values(2 << 1, 2 << 2, 2 << 3, 2 << 4, 2 << 5, 2 << 6, 2 << 7));

INSTANTIATE_TEST_SUITE_P(
    various_Parameters,
    TestPtxtBGV,
    ::testing::Values(BGVParameters(17, 2, 1),
                      BGVParameters(17, 2, 3),
                      BGVParameters(168, 13, 1),
                      BGVParameters(126, 127, 1),
                      BGVParameters(78, 79, 1),
                      BGVParameters(33, 19, 2),
                      // NOTE: This was used because it has 3 good dimensions
                      // BGVParameters(10005, 37, 1),
                      BGVParameters(50, 53, 1)));

} // namespace
