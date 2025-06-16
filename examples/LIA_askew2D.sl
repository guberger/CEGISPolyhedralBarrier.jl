(set-logic LIA)

(synth-inv inv_fun ((x Int) (y Int)))

(define-fun pre_fun ((x Int) (y Int)) Bool
    (and (>= x 10) (>= y 7) (<= (+ x y) 18)))
(define-fun trans_fun ((x Int) (y Int) (x! Int) (y! Int)) Bool
    (and (>= x 0) (>= y 0) (= x! (+ x 1)) (= y! (+ y 2))))
(define-fun post_fun ((x Int) (y Int)) Bool
    (>= (* x 2) y))

(inv-constraint inv_fun pre_fun trans_fun post_fun)

(check-synth)