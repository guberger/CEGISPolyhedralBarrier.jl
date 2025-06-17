(set-logic LIA)

(synth-inv inv_fun ((x Int) (y Int)))

(define-fun pre_fun ((x Int) (y Int)) Bool
    (and (= x 0) (>= y 0) (<= y 50))
)
(define-fun trans_fun ((x Int) (y Int) (x! Int) (y! Int)) Bool
    (or
        (and (<= x 99) (<= x 50) (= x! (+ x 1)) (= y! y))
        (and (<= x 99) (>= x 51) (= x! (+ x 1)) (= y! (+ y 1)))
    )
)
(define-fun post_fun ((x Int) (y Int)) Bool
    (and (<= x 101) (<= y 101))
)

(inv-constraint inv_fun pre_fun trans_fun post_fun)

(check-synth)