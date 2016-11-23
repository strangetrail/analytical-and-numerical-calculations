 ;;(let* ((string "Hello World!")
 ;;       (c-string (cffi:foreign-funcall "strdup" :string string :pointer)))
 ;;  (unwind-protect (write-line (cffi:foreign-string-to-lisp c-string))
 ;;    (cffi:foreign-funcall "free" :pointer c-string :void))
 ;;  (values))
 ;;(cffi:load-foreign-library "libcrypto.so")
 ;;(cffi:defcfun ("MD5" MD5) :void (string :string) (len :int) (ptr :pointer))
 ;;(let ((string-to-convert "The quick brown fox jumped over the lazy dog's back")
 ;;      (ptr (cffi:foreign-alloc :unsigned-char :count 16)))
 ;;  (md5 string-to-convert (length string-to-convert) ptr)
 ;;  (loop for i from 0 below 16 do
 ;;    (format t "~a" (write-to-string (cffi:mem-ref ptr :unsigned-char i) :base 16)))
 ;;  (cffi:foreign-free ptr))
(load "~/common-lisp/asdf/asdf.lisp")
(asdf:load-system :cffi)
(cffi:load-foreign-library "./libquadpack_test.so")
(cffi:defcfun ("launch_quadpack_test" launch-quadpack-test) :int
  (argc :int) (argv :string) (points :pointer))

;;; Alias in case size_t changes.
(defctype size :unsigned-int)
 
;;; To be set as the CURLOPT_WRITEFUNCTION of every easy handle.
(defcallback easy-write size ((ptr :pointer) (size size)
                                                            (nmemb size) (stream :pointer))
(let ((data-size (* size nmemb)))
  (handler-case
    ;; We use the dynamically-bound *easy-write-procedure* to
    ;; call a closure with useful lexical context.
    (progn (funcall (symbol-value '*easy-write-procedure*)
                                            (foreign-string-to-lisp ptr :count data-size))
                          data-size)         ;indicates success
            ;; The WRITEFUNCTION should return something other than the
                    ;; #bytes available to signal an error.
                            (error () (if (zerop data-size) 1 0)))))



(let ((string-to-convert "The quick brown fox jumped over the lazy dog's back")
      (points (cffi:foreign-alloc :double :count 4 :initial-contents #(1.000000e0 2.100000e0 3.100000e0 4.500000e0))))
  (format t "~t~a~%" string-to-convert)
  (format t "~tLISP float length ~a~%" (cffi:foreign-type-size :float))
  (format t "~tLISP double float length ~a~%" (cffi:foreign-type-size :double))
  (launch-quadpack-test (length string-to-convert) string-to-convert points)
  (loop for i from 0 below 4 do
    (format t "~t[~a]:  ~a~%" i (write-to-string (cffi:mem-ref points :double (* i (cffi:foreign-type-size :double))))))
  (cffi:foreign-free points))
