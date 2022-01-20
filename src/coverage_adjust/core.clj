(ns coverage-adjust.core
  (:require [clojure.java.io :as io]
            [clojure.tools.cli :refer [parse-opts]]
            [clojure.data.csv :as csv]
            [clojure.string :as str]
            [clojure.pprint :refer [pprint]])
  (:gen-class))


(def cli-options
  [["-q" "--sequence-qc BASIC_QC"]
   ["-a" "--abundance ABUNDANCE"]
   ["-s" "--species SPECIES"]
   ["-g" "--genome-size GENOME SIZE"]
   ["-h" "--Help"]])


(defn csv-data->maps
  ([csv-data]
   (map zipmap
        (->> (first csv-data) ;; First row is the header
             (map keyword) ;; Drop if you want string keys instead
             repeat)
        (rest csv-data)))
  ([csv-data headers]
   (map zipmap
        (->> headers
             repeat)
        (rest csv-data))))


(defn index-seq-of-maps [k s]
  "Takes a seq of maps (s) and returns a map that is indexed 
   by the values associated with a key (k)

   eg:

   (index-seq-of-maps :id [{:id \"sample-01\" :a 1 :b 2} {:id \"sample-02\" :a 3 :b 4}])

   => {\"sample-01\" {:id \"sample-01\" :a 1 :b 2}
       \"sample-02\" {:id \"sample-02\" :a 3 :b 4}}"
  (into {} (map #(vec [(k %) %])) s))


(defn map-kv [f m]
  "Applies f to all values of map m"
  (reduce-kv #(assoc %1 %2 (f %3)) {} m))


(defn parse-sequence-qc [sequence-qc-path headers]
  (with-open [reader (io/reader sequence-qc-path)]
    (->> (doall (csv-data->maps (csv/read-csv reader) headers))
         (map (fn [x] (update x :estimated-genome-size-bp #(Long/parseLong %))))
         (map (fn [x] (update x :estimated-depth-coverage #(Float/parseFloat %))))
         (map (fn [x] (update x :total-bases #(Long/parseLong %))))
         (map (fn [x] (update x :average-base-quality #(Float/parseFloat %))))
         (map (fn [x] (update x :bases-above-q30-percent #(Float/parseFloat %))))
         (map (fn [x] (update x :gc-content-percent #(Float/parseFloat %))))
         (index-seq-of-maps :library-id))))


(defn parse-abundance [abundance-path]
  (with-open [reader (io/reader abundance-path)]
    (->> (doall (line-seq reader))
         next
         next
         (map #(str/replace % #"\"" ""))
         (map #(str/replace % #"\\" ""))
         (map #(str/split % #","))
         ;; The following step takes this: [["sample-01" "organism-1" "0.90" "organism-2" "0.10"]
         ;;                                 ["sample-02" "organism-1" "0.85" "organism-2" "0.15"]]
         ;; ...and converts to this: {"sample-01" {"organism-1" 0.90 "organism-2" 0.10}
         ;;                           "sample-02" {"organism-1" 0.85 "organism-2" 0.15}}
         (reduce (fn [acc x] (into acc {(first x) (map-kv #(Float/parseFloat %) (apply hash-map (rest x)))})) {}))))


(defn adjust-sample-coverage [sample-seq-qc sample-abundances target-species-name target-species-genome-size]
  (if (contains? sample-abundances target-species-name)
    (let [abundance-percent (get sample-abundances target-species-name)
          adjusted-total-bases (* (/ abundance-percent 100) (:total-bases sample-seq-qc))
          estimated-coverage (Float/parseFloat (format "%.3f" (/ adjusted-total-bases target-species-genome-size)))]
      (assoc sample-seq-qc :target-species-name target-species-name
                           :target-species-abundance abundance-percent
                           :target-species-genome-size-bp target-species-genome-size
                           :estimated-coverage-for-target-genome estimated-coverage))
    (let [estimated-coverage 0.00
          abundance-percent 0.00]
      (assoc sample-seq-qc :target-species-name target-species-name
                           :target-species-abundance abundance-percent
                           :target-species-genome-size-bp target-species-genome-size
                           :estimated-coverage-for-target-genome estimated-coverage))))


(defn -main [& args]
  "Re-estimate depth of coverage, based on number of sequenced bases, 
   target species and genome size"

  ;; Command-line argument parsing
  (def opts (parse-opts args cli-options))

  (def target-species
    (get-in opts [:options :species]))

  (def genome-size
    (Long/parseLong (get-in opts [:options :genome-size])))

  ;; Input file parsing
  (def sequence-qc-headers [:library-id               :estimated-genome-size-bp
                            :estimated-depth-coverage :total-bases
                            :average-base-quality     :bases-above-q30-percent
                            :gc-content-percent])

  (def sequence-qc 
    (parse-sequence-qc (get-in opts [:options :sequence-qc]) sequence-qc-headers))
  
  (def abundances 
    (parse-abundance (get-in opts [:options :abundance])))

  ;; Prepare output by adjusting coverage calculations for each sample
  (def output 
    (let [sample-ids (keys sequence-qc)
          sample-qcs (map #(get sequence-qc %) sample-ids)
          sample-abundances (map #(get abundances %) sample-ids)]
      (map #(adjust-sample-coverage % %2 target-species genome-size) sample-qcs sample-abundances)))

  ;; Print output to stdout
  (let [columns (conj sequence-qc-headers :target-species-name      :target-species-genome-size-bp
                                          :target-species-abundance :estimated-coverage-for-target-genome)
        headers (map #(str/replace (name %) #"-" "_") columns)
        rows (mapv #(mapv % columns) output)]
    (with-open [writer *out*]
      (csv/write-csv writer (cons headers rows)))))


;; Useful forms for development & debugging. Not evaluated during main runtime.
(comment
  (def args [])
  (def opts {:options {:sequence-qc "test_input/sequence_qc.csv"
                       :abundance "test_input/species_abundance.csv"
                       :species "Influenza A virus"
                       :genome-size "13588"}
             :arguments []
             :summary ""
             :errors nil})
  
  (pprint (parse-sequence-qc (get-in opts [:options :sequence-qc])))
  (pprint (parse-abundance (get-in opts [:options :abundance])))
  
  (adjust-sample-coverage 
   (get sequence-qc "")
   (get abundances "")
   "Influenza A virus"
   13588)
  
  )
