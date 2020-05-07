(load "~/common-lisp/asdf/asdf.lisp")

(asdf:load-system :cffi)


(setq *debugger-hook* nil)


;;; It's necessary to use inline functions inside callbacks
;;;  to prevent extra calls:
;;; The name of translated or compiled GNU Maxima function
;;;  should not contain dashes:
(declaim (inline $XY2DINTEGRAND))

;;; Function for building exported C++ code via external SBCL calls
;;;  to /usr/bin/make:
(defun $BUILDCUSTOMINTEGRAND (codefilesname codeprojectname cxxsuffix)
  (let ()
    (progn (sb-ext:run-program
             "/usr/bin/make"
             `("all"
               ,(format nil "codesrcname=~a" codefilesname)
               ,(format nil "codeprojectname=~a" codeprojectname)
               ,(format nil "codesrcsuffix=~a" cxxsuffix))
             :output *standard-output*))))
