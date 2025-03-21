from sqlalchemy import Column, Integer, String, ForeignKey, Float as SAFloat, Table, DateTime, Boolean, func, SmallInteger, event, Enum, UniqueConstraint, ARRAY, select, JSON
from sqlalchemy.orm import relationship, Mapped, mapped_column, declared_attr, composite
from sqlalchemy.dialects.postgresql import UUID
import uuid
from sqlalchemy.orm import declarative_mixin, DeclarativeBase
from typing import Optional, List
from .tools import tablename2plural
# from sqlalchemy_easy_softdelete.mixin import generate_soft_delete_mixin_class
from sqlalchemy.ext.asyncio import AsyncAttrs
from sqlalchemy.dialects.postgresql import JSONB, insert
from sqlalchemy.ext.hybrid import hybrid_property




# class SoftDeleteMixin(generate_soft_delete_mixin_class(
#     # This table will be ignored by the hook
#     # even if the table has the soft-delete column
#     ignored_tables=[])):
#     # type hint for autocomplete IDE support
#     deleted_at: datetime
    
    
class Tablename:
    @declared_attr.directive
    def __tablename__(self) -> Optional[str]:
        st = ''.join(['_' + i.lower() if i.isupper() else i for i in self.__name__]).lstrip('_')
        return tablename2plural(st)


# Necessary for async processus with server_default needed
class EagerDefaultsMixin:
    @declared_attr.directive
    def __mapper_args__(self):
        return {"eager_defaults": True}
    

class Base(AsyncAttrs, DeclarativeBase):
    pass


@declarative_mixin
class OrmHeader(Tablename, EagerDefaultsMixin):
    # __abstract__ = True
    id = Column(UUID(as_uuid=True), primary_key=True, unique=True, default=uuid.uuid4, name='id')
    created_at = Column(DateTime(timezone=True), server_default=func.now(), name='created_at')
    updated_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), name='updated_at')
    deleted_at = Column(DateTime(timezone=True), nullable=True, name='deleted_at')
    created_by = Column(UUID(as_uuid=True), ForeignKey('users.id', ondelete='SET NULL'), nullable=True, name='created_by')
    updated_by = Column(UUID(as_uuid=True), ForeignKey('users.id', ondelete='SET NULL'), nullable=True, name='updated_by')
    deleted_by = Column(UUID(as_uuid=True), ForeignKey('users.id', ondelete='SET NULL'), nullable=True, name='deleted_by')
    archived = Column(Boolean, nullable=True, name='archived')
    archived_at = Column(DateTime(timezone=True), nullable=True, name='archived_at')
    active = Column(Boolean, default=True, name='active')
    locked = Column(Boolean, default=False, name='locked')    
    
    @declared_attr
    def created_users(cls):
        return relationship('User', foreign_keys=[cls.created_by], uselist=True)

    @declared_attr
    def updated_users(cls):
        return relationship('User', foreign_keys=[cls.updated_by], uselist=True)

    @declared_attr
    def deleted_users(cls):
        return relationship('User', foreign_keys=[cls.deleted_by], uselist=True)
    
    
class User(Base, OrmHeader):
    first_name = Column(String(255))
    last_name = Column(String(255))
    email = Column(String(255), unique=True, name='email')
    last_seen = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), name='last_seen')
    realm_roles = Column(ARRAY(String))
    # TODO : Add join_depth in relationship?? --> Eager Loading : lazy + join_depth in relationship : https://docs.sqlalchemy.org/en/20/orm/self_referential.html#configuring-self-referential-eager-loading
    objects: Mapped[List['Object']] = relationship(back_populates='user', primaryjoin="User.id==Object.user_id")
    session_logs: Mapped[List['SessionLog']] = relationship(back_populates='user', primaryjoin="User.id==SessionLog.user_id")


class SessionLog(Base, OrmHeader):
    user_id = Column(UUID(as_uuid=True), ForeignKey('users.id', ondelete='SET NULL'), name='user_id')   
    user: Mapped['User'] = relationship('User', back_populates='session_logs', foreign_keys=[user_id])
    # shiny_id = Column(Integer, nullable=True)
    last_seen = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), name='last_seen')
    # authorization = Column(JSONB, unique=True, name='authorization')
    token = Column(String, unique=True, name='token')
    refresh_token = Column(String, unique=True, name='refresh_token')
    token_expiration_date = Column(DateTime(timezone=True), name='token_expiration_date')
    
    # Unicity Constraint
    __table_args__ = (
        UniqueConstraint('token', name='uq_token'),
    )


class Object(Base, OrmHeader):
    user_id = Column(UUID(as_uuid=True), ForeignKey('users.id', ondelete='SET NULL'), name='user_id')
    user: Mapped['User'] = relationship('User', back_populates='objects', foreign_keys=[user_id])
    name = Column(String(255))
    # score: Mapped['Score'] = relationship(back_populates='float', primaryjoin="Float.id==Score.float_id")
    # f_rules: Mapped[List['FRule']] = relationship(back_populates='float', primaryjoin="Float.id==FRule.float_id")
    # p_rules: Mapped[List['PRule']] = relationship(back_populates='float', primaryjoin="Float.id==PRule.float_id")
    