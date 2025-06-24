(set-logic LIA)

(synth-inv inv_fun ((x Int)))

(define-fun pre_fun ((x Int)) Bool
    (= x 0))
(define-fun trans_fun ((x Int) (x! Int)) Bool
    (and (< x 5) (= x! (+ x 1))))
(define-fun post_fun ((x Int)) Bool
    (<= x 6))

(inv-constraint inv_fun pre_fun trans_fun post_fun)

(check-synth)