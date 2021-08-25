#include <iostream>

#include <helib/helib.h>
#include <helib/EncryptedArray.h>
#include <helib/ArgMap.h>
#include <NTL/BasicThreadPool.h>

// Utility function to read <K,V> CSV data from file
std::vector<std::pair<std::string, int>> read_csv(std::string filename)
{
  std::vector<std::pair<std::string, int>> dataset;
  std::ifstream data_file(filename);

  if (!data_file.is_open())
    throw std::runtime_error(
        "Error: This example failed trying to open the data file: " + filename +
        "\n           Please check this file exists and try again.");

  std::vector<std::string> row;
  std::string line, entry, temp;

  if (data_file.good()) {
    // Read each line of file
    while (std::getline(data_file, line)) {
      row.clear();
      std::stringstream ss(line);
      while (getline(ss, entry, ',')) {
        row.push_back(entry);
      }
      // Add key value pairs to dataset
      std::cout << row[0] << ", " << row[1] << "\n";
      dataset.push_back(std::make_pair(row[0], std::stoi(row[1])));
    }
  }

  data_file.close();
  return dataset;
}

std::vector<std::pair<helib::Ctxt, helib::Ctxt>> setup_database(
    std::vector<std::pair<std::string, int>> db,
    const helib::PubKey& public_key,
    helib::Context& context,
    bool debug = true)
{
  // Convert strings into numerical vectors
  std::cout << "\n---Initializing the encrypted key,value pair database ("
            << db.size() << " entries)...";
  std::cout
      << "\nConverting strings to numeric representation into Ptxt objects ..."
      << std::endl;

  // Generating the Plain text representation of Country DB
  HELIB_NTIMER_START(timer_PtxtCountryDB);
  std::vector<std::pair<helib::Ptxt<helib::BGV>, helib::Ptxt<helib::BGV>>>
      db_ptxt;

  for (const auto& key_value_pair : db) {
    if (debug) {
      std::cout << "\t\tname_addr_pair.first size = "
                << key_value_pair.first.size() << " (" << key_value_pair.first
                << ")"
                << "\tname_addr_pair.second"
                << "(" << key_value_pair.second << ")" << std::endl;
    }

    helib::Ptxt<helib::BGV> country(context);
    std::cout << "\tname size = " << country.size() << std::endl;
    for (long i = 0; i < key_value_pair.first.size(); ++i) {
      country.at(i) = key_value_pair.first[i];
    }

    helib::Ptxt<helib::BGV> capital(context);
    capital.addConstant(key_value_pair.second);

    db_ptxt.emplace_back(std::move(country), std::move(capital));
  }
  HELIB_NTIMER_STOP(timer_PtxtCountryDB);

  // Encrypt the Country DB
  std::cout << "Encrypting the database..." << std::endl;
  HELIB_NTIMER_START(timer_CtxtCountryDB);
  std::vector<std::pair<helib::Ctxt, helib::Ctxt>> encrypted_db;
  for (const auto& key_value_pair : db_ptxt) {
    helib::Ctxt encrypted_key(public_key);
    helib::Ctxt encrypted_value(public_key);
    public_key.Encrypt(encrypted_key, key_value_pair.first);
    public_key.Encrypt(encrypted_value, key_value_pair.second);
    encrypted_db.emplace_back(std::move(encrypted_key),
                              std::move(encrypted_value));
  }
  HELIB_NTIMER_STOP(timer_CtxtCountryDB);

  return encrypted_db;
}

helib::Ctxt length_encrypted_db_key(helib::Ctxt encrypted_key, unsigned long p)
{
  encrypted_key.power(p - 1);
  totalSums(encrypted_key);
  return encrypted_key;
}

int main(int argc, char* argv[])
{
  /************ HElib boiler plate ************/

  // Note: The parameters have been chosen to provide a somewhat
  // faster running time with a non-realistic security level.
  // Do Not use these parameters in real applications.

  // Plaintext prime modulus
  unsigned long p = 131;
  // Cyclotomic polynomial - defines phi(m)
  unsigned long m = 130; // this will give 48 slots
  // Hensel lifting (default = 1)
  unsigned long r = 1;
  // Number of bits of the modulus chain
  unsigned long bits = 1000;
  // Number of columns of Key-Switching matrix (default = 2 or 3)
  unsigned long c = 2;
  // Size of NTL thread pool (default =1)
  unsigned long nthreads = 1;
  // input database file name
  std::string db_filename = "./string_int_dataset.csv";
  // debug output (default no debug output)
  bool debug = true;

  helib::ArgMap amap;
  amap.arg("m", m, "Cyclotomic polynomial ring");
  amap.arg("p", p, "Plaintext prime modulus");
  amap.arg("r", r, "Hensel lifting");
  amap.arg("bits", bits, "# of bits in the modulus chain");
  amap.arg("c", c, "# fo columns of Key-Switching matrix");
  amap.arg("nthreads", nthreads, "Size of NTL thread pool");
  amap.arg("db_filename",
           db_filename,
           "Qualified name for the database filename");
  amap.toggle().arg("-debug", debug, "Toggle debug output", "");
  amap.parse(argc, argv);

  // set NTL Thread pool size
  if (nthreads > 1)
    NTL::SetNumThreads(nthreads);

  std::cout << "\n*********************************************************";
  std::cout << "\n*           Privacy Preserving Search Example           *";
  std::cout << "\n*           =================================           *";
  std::cout << "\n*                                                       *";
  std::cout << "\n* This is a sample program for education purposes only. *";
  std::cout << "\n* It implements a very simple homomorphic encryption    *";
  std::cout << "\n* based db search algorithm for demonstration purposes. *";
  std::cout << "\n*                                                       *";
  std::cout << "\n*********************************************************";
  std::cout << "\n" << std::endl;

  std::cout << "---Initialising HE Environment ... ";
  // Initialize context
  // This object will hold information about the algebra used for this scheme.
  std::cout << "\nInitializing the Context ... ";
  HELIB_NTIMER_START(timer_Context);
  helib::Context context = helib::ContextBuilder<helib::BGV>()
                               .m(m)
                               .p(p)
                               .r(r)
                               .bits(bits)
                               .c(c)
                               .build();
  HELIB_NTIMER_STOP(timer_Context);

  // Secret key management
  std::cout << "\nCreating Secret Key ...";
  HELIB_NTIMER_START(timer_SecKey);
  // Create a secret key associated with the context
  helib::SecKey secret_key = helib::SecKey(context);
  // Generate the secret key
  secret_key.GenSecKey();
  HELIB_NTIMER_STOP(timer_SecKey);

  // Compute key-switching matrices that we need
  HELIB_NTIMER_START(timer_SKM);
  helib::addSome1DMatrices(secret_key);
  HELIB_NTIMER_STOP(timer_SKM);

  // Public key management
  // Set the secret key (upcast: FHESecKey is a subclass of FHEPubKey)
  std::cout << "\nCreating Public Key ...";
  HELIB_NTIMER_START(timer_PubKey);
  const helib::PubKey& public_key = secret_key;
  HELIB_NTIMER_STOP(timer_PubKey);

  // Get the EncryptedArray of the context
  const helib::EncryptedArray& ea = context.getEA();

  // Print the context
  std::cout << std::endl;
  if (debug) {
    context.printout();
  }

  // Print the security level
  // Note: This will be negligible to improve performance time.
  std::cout << "\n***Security Level: " << context.securityLevel()
            << " *** Negligible for this example ***" << std::endl;

  // Get the number of slot (phi(m))
  long nslots = ea.size();
  std::cout << "\nNumber of slots: " << nslots << std::endl;

  /************ Read in the database ************/
  std::vector<std::pair<std::string, int>> country_db;
  try {
    country_db = read_csv(db_filename);
  } catch (std::runtime_error& e) {
    std::cerr << "\n" << e.what() << std::endl;
    exit(1);
  }

  std::vector<std::pair<helib::Ctxt, helib::Ctxt>> encrypted_country_db =
      setup_database(country_db, public_key, context);

  // Print DB Creation Timers
  if (debug) {
    helib::printNamedTimer(std::cout << std::endl, "timer_Context");
    helib::printNamedTimer(std::cout, "timer_Chain");
    helib::printNamedTimer(std::cout, "timer_SecKey");
    helib::printNamedTimer(std::cout, "timer_SKM");
    helib::printNamedTimer(std::cout, "timer_PubKey");
    helib::printNamedTimer(std::cout, "timer_PtxtCountryDB");
    helib::printNamedTimer(std::cout, "timer_CtxtCountryDB");
  }

  std::cout << "\nInitialization Completed - Ready for Queries" << std::endl;
  std::cout << "--------------------------------------------" << std::endl;

  /** Create the query **/

  // Read in query from the command line
  std::string query_string;
  std::cout << "\nPlease enter a valid name: ";
  // std::cin >> query_string;
  std::getline(std::cin, query_string);
  std::cout << "Looking for the value of name: " << query_string << std::endl;
  std::cout << "This may take few minutes ... " << std::endl;

  HELIB_NTIMER_START(timer_TotalQuery);

  HELIB_NTIMER_START(timer_EncryptQuery);

  // Convert query to a numerical vector
  helib::Ptxt<helib::BGV> query_ptxt(context);
  for (long i = 0; i < query_string.size(); ++i) {
    query_ptxt[i] = query_string[i];
  }

  // Encrypt the query
  helib::Ctxt query(public_key);
  public_key.Encrypt(query, query_ptxt);
  HELIB_NTIMER_STOP(timer_EncryptQuery);

  /************ Perform the database search ************/

  HELIB_NTIMER_START(timer_QuerySearch);
  std::vector<helib::Ctxt> mask;
  mask.reserve(country_db.size());
  for (const auto& encrypted_pair : encrypted_country_db) {
    helib::Ctxt mask_entry = encrypted_pair.first; // Copy of database key

    helib::Ctxt mask_entry_length = length_encrypted_db_key(mask_entry, p);

    mask_entry -= query;                 // Calculate the difference
    mask_entry.power(p - 1);             // Fermat's little theorem
    mask_entry.negate();                 // Negate the ciphertext
    mask_entry.addConstant(NTL::ZZX(1)); // 1 - mask = 0 or 1

    // Create a vector of copies of the mask
    std::vector<helib::Ctxt> rotated_masks(ea.size(), mask_entry);
    for (int i = 1; i < rotated_masks.size(); i++) {
      ea.rotate(rotated_masks[i], i); // Rotate each of the masks
    }

    totalProduct(mask_entry, rotated_masks);      // Multiply each of the masks
    mask_entry.multiplyBy(encrypted_pair.second); // multiply mask with values
    mask.push_back(mask_entry);
  }

  // Aggregate the results into a single ciphertext
  // Note: This code is for educational purposes and thus we try to refrain
  // from using the STL and do not use std::accumulate
  helib::Ctxt value = mask[0];
  for (int i = 1; i < mask.size(); i++) {
    value += mask[i];
  }

  HELIB_NTIMER_STOP(timer_QuerySearch);

  /************ Decrypt and print result ************/

  HELIB_NTIMER_START(timer_DecryptQueryResult);
  helib::Ptxt<helib::BGV> plaintext_result(context);
  secret_key.Decrypt(plaintext_result, value);
  HELIB_NTIMER_STOP(timer_DecryptQueryResult);

  // extract first element as the result
  std::string query_result =
      std::to_string(static_cast<long>(plaintext_result[0]));
  HELIB_NTIMER_STOP(timer_TotalQuery);

  if (query_result == "0") {
    query_result = "Name not in the database."
                   "\n*** Please make sure to enter a valid name"
                   "\n*** with the first letter in upper case.";
  }
  std::cout << "\nQuery result: " << query_result << std::endl;

  // Print DB Query Timers
  if (debug) {
    helib::printNamedTimer(std::cout << std::endl, "timer_EncryptQuery");
    helib::printNamedTimer(std::cout, "timer_QuerySearch");
    helib::printNamedTimer(std::cout, "timer_DecryptQueryResult");
    std::cout << std::endl;
  }

  helib::printNamedTimer(std::cout, "timer_TotalQuery");

  return 0;
}
