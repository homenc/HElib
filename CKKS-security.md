# Security of Approximate-Numbers Homomorphic Encryption

Security of standard (not homomorphic) encryption scheme is fairly well understood nowadays, and it usually comes in two flavors: An encryption scheme may only resists passive attackers (that can view encrypted data but not manipulate it), or it can also resist active attackers that can manipulate the encrypted data and then observe how it is decrypted.
The crypto lingo for the weaker (passive) notion is [CPA-security](https://en.wikipedia.org/wiki/Chosen-plaintext_attack), while the stronger is called [CCA-security](https://en.wikipedia.org/wiki/Chosen-ciphertext_attack). (These acronyms stand for chosen-plaintext-attacks and chosen-ciphertext-attacks, respectively.)
It is well understood that to protect data at rest or in transit, one must use CCA-secure encryption. CPA-security may suffice when the encryption is used as one component in a larger system, and protection against active attackers is provided by other components in that system.

For homomorphic encryption (HE) schemes, it is well known that they inherently cannot be CCA secure. One consequence of their non-CCA-security is that the system must ensure that decryption is never applied to invalid ciphertexts.<sup>[1](#validCtxt)</sup> Indeed, with contemporary HE schemes, allowing the attacker to submit invalid ciphertexts to be decrypted would typically result in exposure of the secret key. It is thus a commonly accepted practice for HE schemes to make do with CPA security, and rely on the system around the HE to provide the extra protection that may be needed.

Recently, Li and Micciancio [observed](https://eprint.iacr.org/2020/1533) that for *approximate-number* HE schemes, the common notion of CPA security may not be enough even against passive attackers. The crucial difference is that since the scheme itself adds some error, the attacker could learn something from the decryption result *even if it knows what the decrypted result was supposed to be*. Specifically, it learns the error, which may leak information about the secret key. Li and Micciancio described simple attacks on the CKKS approximate-numbers scheme that expose the secret key after seeing a small number of decryption results (sometimes as few as a single decryption).

## Mitigating the Li-Micciancio attacks

In response to this observation, HElib now includes countermeasures to mitigate this risk. In particular, the HElib decryption routine for CKKS was modified to add some key-independent noise, masking the key-dependent noise and making the Li-Micciancio attacks harder to mount. This extra noise, however, reduces the accuracy of the decrypted result. By default, the magnitude of this key-independent noise is chosen equal to the noise-bound that HElib keeps for the ciphertext to be decrypted, but that bound may be over-conservative and hence the added noise could sometimes cause a significant decrease in accuracy. HElib therefore provides the application the means to specify the required accuracy on decryption, and will then add the largest amount of noise subject to this accuracy constraint.

Applications that use CKKS need to walk a fine line, balancing security, accuracy, and performance. Below we describe a framework for developing CCKS-based applications while taking these aspect into account. Roughly, we recommend that the application try to get a tight estimate for the error, and then use that estimate instead of the noise bound that HElib computes.

1. Start by identifying the (class of) processes that the application needs to apply to the data;
2. Identify the *required precision* from the processed results;
3. Run these processes in HElib on test data to determine the parameters to be used, roughly as follows:
   1. Find some parameter setting where `context.securityLevel()` returns high enough security and the library does not report decryption errors;
   2. Experimentally find a fairly tight bound on the noise magnitude of the decrypted results, when called with the accuracy parameter from Step 2. This can be done by running HElib processing many times and comparing the decrypted results to what you get when applying the same processing to plaintext data;
   3. Compare the noise bounds that you get from Step 3.2 above with the estimated error that HElib computes just prior to decryption, which can be accessed via `ctxt.errorBound()`. If the experimental noise is significantly smaller than what `ctxt.errorBound()` returns then use `ctxt.bumpNoiseBound()` prior to decryption to force HElib to use the tighter error estimate.
   4. Repeat Steps 3.1-3.3, increasing the parameters until decryption no longer emits a warning about the added noise being too small.

We now elaborate more on each of these steps.

### Steps 1-2, Processes and Required Precision

As with any application of HE, the starting point is ensuring adequate functionality. The developers must therefore determine the operations to be computed and the required *output precision*. Some applications only need to compute a single function while others need to handle a wider class of functions, but at the very least the developers should determine the maximum depth of supported circuits and the maximum precision required.

In HElib, the application gets to specify the *input precision* for encrypting plaintext data, and each level in the circuit reduces the precision by one (or upto 1.5) bits, namely the error roughly doubles for each level. For example using the helper class `PtxtArray` one could call `ptxt.encrypt(ctxt, ptxtMagnitude, precision)`,<sup>[2](#optPrecision)</sup> which would result in encryption of the given input plaintext with error bounded by 2<sup>-precision</sup>.
After evaluating a depth-three circuit (say), the error will grow to a little more than 2<sup>3-precision</sup>. An application requiring *p* bits of output precision after the evaluation of a depth-*d* circuit should then encrypt with input precision somewhere between *p+d* and *p+3d/2*.

### Step 3, Finding Matching Parameters

After determining the required depth *d* and output precision *p*, the application developers need to find some parameters that yield this output precision while ensuring security. This typically involves an iterative trial-and-error procedure, as follows:

1. Begin with some rough estimate of the parameters: The total bitlength of the modulus needed for a depth-*d* circuit is typically something like *22(d+1)*, and the ring dimension needed to get 128-bit security with this modulus is around *750(d+1)* (rounded up to a power of two, since the current CKKS support in HElib was only tested with power-of-two rings).

   Denote by *m* the next power of two larger than *750(d+1)*, then these initial parameters can be set in HElib by using the `ContextBuilder` class, setting
   ```
     Context context =
         ContextBuilder<CKKS>().m(m).bits(16*(d+1)).precision(p+d).c(3).build();
   ```
   Note that the number of bits specified above is only *16(d+1)* rather than *22(d+1)*, and there is a mysterious additional parameter *c=3* satisfying *(c+1)/c ~ 22/16*. The quantity *bits=16(d+1)* is an estimate for what's needed for a depth-*d* circuit, and the total modulus bitsize in HElib is obtained roughly as *bits&middot;(c+1)/c*. (See discussion of *ciphertext vs. special primes* in the [HElib document](https://eprint.iacr.org/2020/1481), section 5.1.) One can get slightly smaller total bitsize for the same functionality by increasing the *c* parameter from its default value of *c=3*, but the keysize and running time grow almost linearly with *c*.

   After setting these parameters, run HElib processing and check if the library complains about possible decryption errors, and use these experiments to increase or decrease the *bits* and *m* parameters. You should also use the interface `context.securityLevel()` to check that this combination of parameters *m,bits,precision* and *c* yields an acceptable level of security.

2. Check that these parameters indeed yield the required output precision. Namely run several experiments of HElib processing, decrypting the result with the optional precision argument (e.g., using the interface `ptxt.rawDecrypt(ctxt, secretKey)` of the `ptxtArray` class), and comparing the outcome with the result of plaintext processing. (At this point you should ignore any warning from HEib about the added noise being too small.) Revisit the parameters from Step 1 above until both security and output precision are acceptable.

3. Next, check the error bound that HElib reports via `ctxt.errorBound()` just prior to decryption, to the actual error that was observed in the previous step (which presumably is no more than *2<sup>-p</sup>*). If the HElib estimate is significantly larger than the actual noise, then use `ctxt.bumpNoiseBound()` prior to decryption to force HElib to use the "right estimate".

4. At this point, you need to switch to using `ptxt.decrypt(ctxt, secretKey, p)` and check that HElib no longer emits a warning message on decryption about the added noise being too small. If you still see that warning, then you need to revisit the parameter setting, increasing the *bits* parameter (and other parameters as needed to maintain security and precision), until this warning is no longer displayed.

As a final comment, we remark that this procedure is only adequate for cases where the threat model include the adversary only getting access to a handful of decryption results. If the total amount of information available to the adversary about decryption queries is much larger, then the application needs to increase the *bits* parameter from above. In particular, to withstand an attack that has access to the results from decryption of *D* ciphertexts, the *bits* parameter should be increased roughly by a factor of *<math><msqrt><mi>D</mi></msqrt></math>*, and the input-precision parameter should similarly be increased to *<math><msqrt><mi>D</mi></msqrt>(p+d)</math>*.

---------------------------------------
<a name="validCtxt"><sup>1</sup></a>
In HElib, a ciphertext is *valid* if it was produced using the library's encryption and evaluation routines and does not contain too much noise. In particular, decryption using should not emit a warning about "decrypting with too much noise".

<a name="optPrecision"><sup>2</sup></a>
The `ptxtMagnitude`, and `precision` arguments are optional, with defaults that depend on the context and the actual plaintext to be encrypted. Importantly, the precision (if specified) indicates an absolute number, *not* a fraction of the plaintext magnitude. For example calling with `precision=3` yields an error of magnitude 1/8, whether it is called with `ptxtMagnitude=16`, and `ptxtMagnitude=1`. Finally note that the `precision` arguments can also be negative, indicating error magnitude larger than one.

