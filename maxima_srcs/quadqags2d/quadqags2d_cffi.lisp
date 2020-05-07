(cffi:load-foreign-library "quadqags2d/libquadqags2d.so")

;;; Definition of wrapped call of resource management functions:
(cffi:defcfun ("allocGIW" alloc-giw) :void)
(cffi:defcfun ("freeGIW" free-giw) :void)

;;; Definition of wrapped call of the main library function:
(cffi:defcfun ("quad_qags_2d" quad-qags-2d) :double
  (l_lower_x :double) (l_lower_y :double)
  (l_upper_x :double) (l_upper_y :double)
  (params :pointer)
  (integrand :pointer))

;;; Alias in case `size_t' changes:
(cffi:defctype size :unsigned-int)

;;; Callback function:
(cffi:defcallback integrand-2d :double ((x :double) (y :double)
                                        (rawparams :pointer))
  (let ((xx (cffi:convert-from-foreign x :double))
        (yy (cffi:convert-from-foreign y :double))
        (theomega (cffi:convert-from-foreign (cffi:mem-ref
                                               rawparams
                                               :double
                                               0)
                    :double))
        (zz (cffi:convert-from-foreign (cffi:mem-ref
                                         rawparams
                                         :double
                                         (cffi:foreign-type-size :double))
              :double))
        (tt (cffi:convert-from-foreign (cffi:mem-ref
                                         rawparams
                                         :double
                                         (* 2 (cffi:foreign-type-size
                                                :double)))
              :double))
        (thephi0 (cffi:convert-from-foreign (cffi:mem-ref
                                              rawparams
                                              :double
                                              (* 3 (cffi:foreign-type-size
                                                     :double)))
                   :double))
        (thealphaeuler (cffi:convert-from-foreign (cffi:mem-ref
                                                    rawparams
                                                    :double
                                                    (* 4
                                                      (cffi:foreign-type-size
                                                        :double)))
                         :double)))
    (progn (cffi:convert-to-foreign ($XY2DINTEGRAND xx yy theomega zz tt
                                      thephi0 thealphaeuler)
            :double))))

;;; GSL resource management routines:
(defmfun $ALLOCGIW () (let () (progn (alloc-giw))))
(defmfun $FREEGIW () (let () (progn (free-giw))))

;;; 2D QAGS C GSL alternative to GNU Maxima internal QUADPACK implementation:
(defmfun $QUADQAGS2D (lx ly ux uy parameters)
  (let ((lower-boundary-x (cffi:convert-to-foreign lx :double))
        (lower-boundary-y (cffi:convert-to-foreign ly :double))
        (upper-boundary-x (cffi:convert-to-foreign ux :double))
        (upper-boundary-y (cffi:convert-to-foreign uy :double))
        (params (cffi:foreign-alloc :double
                  :count (array-dimension parameters 0)
                  :initial-contents parameters))
        (qags-result 0.0e0))
    (setq qags-result (progn (cffi:convert-from-foreign
                              (quad-qags-2d
                                lower-boundary-x
                                lower-boundary-y
                                upper-boundary-x
                                upper-boundary-y
                                params
                                (cffi:callback integrand-2d))
                              :double)))
    (cffi:foreign-free params)
    (progn qags-result)))
