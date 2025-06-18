(set-logic LIA)

(synth-inv InvF ((T Int) (t Int) (p Bool) (q Bool)))

(define-fun PreF ((T Int) (t Int) (p Bool) (q Bool)) Bool
    (and (>= T 16) (<= T 24) (= t 0) (not q))
)

(define-fun TransF ((T Int) (t Int) (p Bool) (q Bool)
                    (T! Int) (t! Int) (p! Bool) (q! Bool)) Bool
    (and
        (or
            (and p (>= T 19) (= T! (- T 1)) p!)
            (and p (<= T 18) (= T! (+ T 1)) (not p!))
            (and (not p) (<= T 21) (= T! (+ T 1)) (not p!))
            (and (not p) (>= T 22) (= T! (- T 1)) p!)
        )
        (or
            (and (not q) (<= t 48) (= t! (+ t 1)) (not q!))
            (and (>= t 49) (= t! 50) q!)
        )
    )
)

(define-fun PostF ((T Int) (t Int) (p Bool) (q Bool)) Bool
    (and
        (>= t -1) (<= t 51)
        (> T 0) (< T 30)
        (or (not q) (and (>= T 17) (<= T 23)))
    )
)

(inv-constraint InvF PreF TransF PostF)

(check-synth)