# -*- coding: utf-8 -*-
# Generated by the protocol buffer compiler.  DO NOT EDIT!
# source: mixing_state.proto
"""Generated protocol buffer code."""
from google.protobuf.internal import builder as _builder
from google.protobuf import descriptor as _descriptor
from google.protobuf import descriptor_pool as _descriptor_pool
from google.protobuf import symbol_database as _symbol_database
# @@protoc_insertion_point(imports)

_sym_db = _symbol_database.Default()


import matrix_pb2 as matrix__pb2


DESCRIPTOR = _descriptor_pool.Default().AddSerializedFile(b'\n\x12mixing_state.proto\x12\x08\x62\x61yesmix\x1a\x0cmatrix.proto\"\x1c\n\x07\x44PState\x12\x11\n\ttotalmass\x18\x01 \x01(\x01\"-\n\x07PYState\x12\x10\n\x08strength\x18\x01 \x01(\x01\x12\x10\n\x08\x64iscount\x18\x02 \x01(\x01\"9\n\nLogSBState\x12+\n\x11regression_coeffs\x18\x01 \x01(\x0b\x32\x10.bayesmix.Matrix\"V\n\x0cTruncSBState\x12 \n\x06sticks\x18\x01 \x01(\x0b\x32\x10.bayesmix.Vector\x12$\n\nlogweights\x18\x02 \x01(\x0b\x32\x10.bayesmix.Vector\")\n\x08MFMState\x12\x0e\n\x06lambda\x18\x01 \x01(\x01\x12\r\n\x05gamma\x18\x02 \x01(\x01\"\x98\x02\n\x0bMixingState\x12%\n\x08\x64p_state\x18\x01 \x01(\x0b\x32\x11.bayesmix.DPStateH\x00\x12%\n\x08py_state\x18\x02 \x01(\x0b\x32\x11.bayesmix.PYStateH\x00\x12,\n\x0clog_sb_state\x18\x03 \x01(\x0b\x32\x14.bayesmix.LogSBStateH\x00\x12\x30\n\x0etrunc_sb_state\x18\x04 \x01(\x0b\x32\x16.bayesmix.TruncSBStateH\x00\x12\'\n\tmfm_state\x18\x05 \x01(\x0b\x32\x12.bayesmix.MFMStateH\x00\x12)\n\rgeneral_state\x18\x06 \x01(\x0b\x32\x10.bayesmix.VectorH\x00\x42\x07\n\x05stateb\x06proto3')

_builder.BuildMessageAndEnumDescriptors(DESCRIPTOR, globals())
_builder.BuildTopDescriptorsAndMessages(DESCRIPTOR, 'mixing_state_pb2', globals())
if _descriptor._USE_C_DESCRIPTORS == False:

  DESCRIPTOR._options = None
  _DPSTATE._serialized_start=46
  _DPSTATE._serialized_end=74
  _PYSTATE._serialized_start=76
  _PYSTATE._serialized_end=121
  _LOGSBSTATE._serialized_start=123
  _LOGSBSTATE._serialized_end=180
  _TRUNCSBSTATE._serialized_start=182
  _TRUNCSBSTATE._serialized_end=268
  _MFMSTATE._serialized_start=270
  _MFMSTATE._serialized_end=311
  _MIXINGSTATE._serialized_start=314
  _MIXINGSTATE._serialized_end=594
# @@protoc_insertion_point(module_scope)
