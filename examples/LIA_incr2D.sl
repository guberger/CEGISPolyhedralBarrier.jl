(set-logic LIA)

(synth-inv inv_fun ((x Int) (y Int)))

(define-fun pre_fun ((x Int) (y Int)) Bool
    (and (>= x -5) (<= x 5) (>= y -5) (<= y 5)))
(define-fun trans_fun ((x Int) (y Int) (x! Int) (y! Int)) Bool
    (and (= x! (+ x 1)) (= y! (+ y 2))))
(define-fun post_fun ((x Int) (y Int)) Bool
    (>= (* x 2) (- y 15)))

(inv-constraint inv_fun pre_fun trans_fun post_fun)

(check-synth)