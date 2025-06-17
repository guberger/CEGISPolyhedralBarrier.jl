(set-logic LIA)

(synth-inv inv_fun ((x Int) (y Int) (t Int) (q Bool)))

(define-fun pre_fun ((x Int) (y Int) (t Int) (q Bool)) Bool
    (and (>= x -10) (<= x 10) (>= y -10) (<= y 10) (= t 0) (not q))
)
(define-fun trans_fun ((x Int) (y Int) (t Int) (q Bool)
                       (x! Int) (y! Int) (t! Int) (q! Bool)) Bool
    (and
        (or
            (and (>= x (+ y 1)) (= x! (- x 1)) (= y! y))
            (and (>= y (+ x 1)) (= x! x) (= y! (- y 1)))
            (and (= x y) (= x! x) (= y! y))
        )
        (or
            (and (<= t 18) (not q) (= t! (+ t 1)) (not q!))
            (and (>= t 19) (= t! 20) q!)
        )
    )
)
(define-fun post_fun ((x Int) (y Int) (t Int) (q Bool)) Bool
    (or
        (and (not q) (>= x -11) (<= x 11) (>= y -11) (<= y 11))
        (and q (>= x -11) (<= x 11) (>= y -11) (<= y 11) (<= x y))
    )
)

(inv-constraint inv_fun pre_fun trans_fun post_fun)

(check-synth)