#ifndef CELLIDDECODER_H
#define CELLIDDECODER_H

// Standard
#include <string>
#include <iostream>
#include <map>

// ACTSTracking
#include "BitField64.hxx"

/**
 * @class CellIDDecoder
 * @brief Decoder for cell IDs using a bitfield.
 * This is adapted from LCIO's CellIDDecoder to work with
 * edm4hep. It is very similar
 * @author Samuel Ferraro
 */
class CellIDDecoder {
    public:
	CellIDDecoder() : _field("") {} ///<Default constructor
	/**
         * @brief Constructor for CellIDDecoder.
         * @param initString Initialization string describing the fields.
         */
	CellIDDecoder(const std::string& initString) : _field(initString) {}
        /**
         * @brief Template constructor for CellIDDecoder.
         * @param coll Collection (unused in this implementation).
         * @param initString Initialization string describing the fields.
         */
	template<typename T>
        CellIDDecoder(const T& coll, const std::string& initString) : _field(initString) {}
	/**
         * @brief Accesses a field by name.
         * @param name The name of the field.
         * @return Const reference to the BitFieldValue.
         */
        const BitFieldValue& operator()(const std::string& name) const {
            return _field[name];
        }
	/**
         * @brief Gets a description of the fields in the bitfield.
         * @return A string describing the fields.
         */
        std::string fieldDescription() const {
            return _field.fieldDescription();
        }
	/**
         * @brief Sets the value of the bitfield.
         * @param val The value to set.
         */
        void setValue(unsigned val) {
            _field.setValue(val);
        }
	/**
         * @brief Resets the bitfield to zero.
         */
        void reset() {
            _field.reset();
        }

    private:
        mutable BitField64 _field; ///< Bitfield used for decoding.
    };

#endif // CELLIDDECODER_H
